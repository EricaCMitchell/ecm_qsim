import math

import pytest

from ecm_qsim import brickwork, sim

needs_sim = pytest.mark.skipif(
    not sim.default_sim_path().exists(), reason="sim binary not built"
)


def lines(text: str) -> list[str]:
    """The circuit lines after the header and the QUBITS line."""
    return text.splitlines()[2:]


def steps(text: str) -> list[list[str]]:
    """The lines of each time step, split at TICK."""
    result = [[]]
    for line in lines(text):
        if line == "TICK":
            result.append([])
        else:
            result[-1].append(line)
    # The final TICK closes the last step and adds no step of its own.
    assert result[-1] == []
    return result[:-1]


def test_header_records_parameters():
    text = brickwork.brickwork_circuit(4, 3, seed=11)
    assert text.splitlines()[:2] == ["# brickwork n=4 depth=3 seed=11", "QUBITS 4"]


def test_circuit_has_two_steps_per_layer():
    text = brickwork.brickwork_circuit(5, 7, seed=1)
    assert text.endswith("TICK\n")
    assert len(steps(text)) == 14


def test_every_rotation_step_has_ry_on_each_qubit():
    s = steps(brickwork.brickwork_circuit(5, 4, seed=1))
    for step in s[0::2]:
        assert [line.split()[0] for line in step] == ["RY"] * 5
        assert [int(line.split()[2]) for line in step] == [0, 1, 2, 3, 4]
        assert all(0 <= float(line.split()[1]) < 2 * math.pi for line in step)


def test_first_layer_starts_with_pair_0_1():
    s = steps(brickwork.brickwork_circuit(6, 2, seed=1))
    assert s[1] == ["CNOT 0 1", "CNOT 2 3", "CNOT 4 5"]
    assert s[3] == ["CNOT 1 2", "CNOT 3 4"]


def test_odd_n_leaves_last_qubit_unpaired():
    s = steps(brickwork.brickwork_circuit(5, 2, seed=1))
    assert s[1] == ["CNOT 0 1", "CNOT 2 3"]
    assert s[3] == ["CNOT 1 2", "CNOT 3 4"]


def test_two_qubits_keep_empty_cnot_steps():
    s = steps(brickwork.brickwork_circuit(2, 2, seed=1))
    assert len(s) == 4
    assert s[1] == ["CNOT 0 1"] and s[3] == []


def test_same_seed_gives_same_circuit():
    assert brickwork.brickwork_circuit(4, 5, seed=3) == brickwork.brickwork_circuit(
        4, 5, seed=3
    )
    assert lines(brickwork.brickwork_circuit(4, 5, seed=3)) != lines(
        brickwork.brickwork_circuit(4, 5, seed=4)
    )


def test_angles_round_trip_exactly():
    text = brickwork.brickwork_circuit(4, 5, seed=2)
    assert "np." not in text
    angles = [line.split()[1] for line in lines(text) if line.startswith("RY")]
    assert len(angles) == 20
    assert all(repr(float(t)) == t for t in angles)


@pytest.mark.parametrize("n, depth", [(1, 20), (4, 0)])
def test_rejects_too_few_qubits_or_layers(n, depth):
    with pytest.raises(ValueError):
        brickwork.brickwork_circuit(n, depth, seed=1)


@needs_sim
def test_circuit_runs_in_sim():
    counts = sim.run(brickwork.brickwork_circuit(4, 20, seed=1), shots=10, seed=1)
    assert sum(counts.values()) == 10
    assert all(len(bits) == 4 for bits in counts)
