import numpy as np
import pytest

from ecm_qsim import idle, sim


def test_idle_circuit_has_x_step_then_m_empty_steps():
    assert idle.idle_circuit(0) == "QUBITS 1\nX 0\nTICK\n"
    assert idle.idle_circuit(3) == "QUBITS 1\nX 0\nTICK\nTICK\nTICK\nTICK\n"
    with pytest.raises(ValueError):
        idle.idle_circuit(-1)


def test_analytic_p_without_idle_steps_is_one_minus_readout():
    assert idle.analytic_p(0, 0.3, 0.02) == pytest.approx(0.98)


def test_analytic_p_after_one_step_flips_with_two_thirds_p_idle():
    # X and Y flip the bit; Z does not.
    assert idle.analytic_p(1, 0.3, 0.0) == pytest.approx(0.8)


def test_analytic_p_decays_to_one_half():
    assert idle.analytic_p(10_000, 0.01, 0.05) == pytest.approx(0.5)


def test_fit_recovers_p_idle_from_analytic_p():
    m = np.arange(0, 200, 10)
    p = idle.analytic_p(m, 0.012, 0.01)
    fit = idle.fit_p_idle(m, p, 10_000, p_readout=0.01)
    assert fit.p_idle == pytest.approx(0.012, rel=1e-6)
    assert 0 < fit.stderr < 0.012 / 50


def test_fit_stderr_matches_spread_of_repeated_fits():
    rng = np.random.default_rng(1)
    m = np.arange(0, 150, 15)
    shots = 2000
    p = idle.analytic_p(m, 0.02, 0.01)
    fits = [
        idle.fit_p_idle(m, rng.binomial(shots, p) / shots, shots, p_readout=0.01)
        for _ in range(300)
    ]
    estimates = np.array([fit.p_idle for fit in fits])
    assert estimates.mean() == pytest.approx(0.02, rel=0.02)
    assert estimates.std() == pytest.approx(fits[0].stderr, rel=0.15)


@pytest.mark.parametrize("p_idle", [1e-4, 0.7])
def test_fit_recovers_p_idle_across_range(p_idle):
    m = idle.default_m_values(p_idle)
    p = idle.analytic_p(m, p_idle, 0.01)
    fit = idle.fit_p_idle(m, p, 10_000, p_readout=0.01)
    assert fit.p_idle == pytest.approx(p_idle, rel=1e-4)


def test_fit_handles_zero_readout_with_every_shot_returning_one():
    m = np.array([0, 50, 100])
    p = np.array([1.0, 0.76, 0.62])
    fit = idle.fit_p_idle(m, p, 1000, p_readout=0.0)
    assert 0 < fit.p_idle < 0.05
    assert 0 < fit.stderr < np.inf


def test_fit_at_zero_p_idle_returns_zero():
    m = np.array([0, 50, 100])
    p = np.full(3, 0.99)
    fit = idle.fit_p_idle(m, p, 10_000, p_readout=0.01)
    assert fit.p_idle == pytest.approx(0, abs=1e-6)
    assert 0 < fit.stderr < np.inf


def test_fit_without_idle_steps_raises():
    with pytest.raises(ValueError):
        idle.fit_p_idle(np.array([0, 0]), np.array([0.99, 0.99]), 1000, p_readout=0.01)


def test_binomial_stderr_is_never_zero():
    stderr = idle.binomial_stderr(np.array([0.0, 0.5, 1.0]), 100)
    assert np.all(stderr > 0)
    assert stderr[1] == pytest.approx(0.05)


def test_default_m_values_span_three_decay_lengths():
    m = idle.default_m_values(0.01)
    assert m[0] == 0
    assert m[-1] == 223  # 3 decay lengths = 223.497 steps
    assert len(m) <= 21
    assert np.all(np.diff(m) > 0)
    for p_idle in (0.0, 0.75):
        with pytest.raises(ValueError):
            idle.default_m_values(p_idle)


@pytest.mark.skipif(not sim.default_sim_path().exists(), reason="sim binary not built")
def test_run_matches_model():
    m = np.array([0, 20, 60])
    data = idle.run(m, p_idle=0.02, p_readout=0.02, shots=20_000, seed=7)
    assert list(data.m) == [0, 20, 60]
    assert list(data.shots) == [20_000] * 3
    # Five standard errors, so a correct simulator fails very rarely.
    expected = idle.analytic_p(m, 0.02, 0.02)
    assert np.all(np.abs(data.p - expected) < 5 * data.p_stderr)
