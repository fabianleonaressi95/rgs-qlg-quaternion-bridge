import numpy as np
import pytest

from src.q6_complex import build_q6_complex
from src.q6_geometric_hodge import (
    audit_cubical_whitney_hodge,
    build_cubical_whitney_hodge,
    build_unit_cube_complex,
    cubical_whitney_mass_matrices,
)


def test_full_six_cube_cell_counts_and_betti_numbers():
    report = audit_cubical_whitney_hodge()
    assert report["cell_counts_by_degree"] == [64, 192, 240, 160, 60, 12, 1]
    assert report["boundary_ranks"] == [63, 129, 91, 35, 7, 1]
    assert report["betti_numbers"] == [1, 0, 0, 0, 0, 0, 0]
    assert report["harmonic_dimensions"] == [1, 0, 0, 0, 0, 0, 0]
    assert report["pass"]



def test_full_cube_one_skeleton_matches_canonical_q6_edge_orientation():
    q6 = build_q6_complex()
    cube = build_unit_cube_complex(6)
    edge_index = {edge: i for i, edge in enumerate(q6.edges)}
    for column, cell in enumerate(cube.cells[1]):
        axis = cell.index(-1)
        lower_bits = [0 if value == -1 else value for value in cell]
        upper_bits = lower_bits.copy()
        upper_bits[axis] = 1
        lower = int("".join(map(str, lower_bits)), 2)
        upper = int("".join(map(str, upper_bits)), 2)
        old_column = edge_index[(lower, upper)]
        assert np.array_equal(cube.boundary[1][:, column], q6.boundary_1[:, old_column])

def test_boundary_squared_vanishes_exactly_in_all_degrees():
    complex_ = build_unit_cube_complex(6)
    for k in range(2, 7):
        residual = complex_.boundary[k - 1].astype(np.int64) @ complex_.boundary[k].astype(np.int64)
        assert np.count_nonzero(residual) == 0


def test_whitney_mass_matrices_are_symmetric_positive_definite():
    masses = cubical_whitney_mass_matrices(build_unit_cube_complex(6))
    for mass in masses:
        assert np.allclose(mass, mass.T, atol=1e-13)
        assert np.linalg.eigvalsh(mass)[0] > 0


def test_full_cube_fills_the_111_harmonic_two_cochains_of_the_two_skeleton():
    # The 2-skeleton has b2=111, but the full cube is contractible.
    q6 = build_q6_complex()
    rank_b1 = np.linalg.matrix_rank(q6.boundary_1.astype(float))
    rank_b2 = np.linalg.matrix_rank(q6.boundary_2.astype(float))
    assert 240 - rank_b2 == 111
    full = audit_cubical_whitney_hodge()
    assert full["betti_numbers"][2] == 0


@pytest.mark.parametrize("dimension", [1, 2, 3, 4])
def test_dimension_generic_small_cubes(dimension):
    report = audit_cubical_whitney_hodge(dimension=dimension)
    assert report["betti_numbers"] == [1] + [0] * dimension
    assert report["pass"]


def test_invalid_dimension_rejected():
    with pytest.raises(ValueError):
        build_unit_cube_complex(0)
