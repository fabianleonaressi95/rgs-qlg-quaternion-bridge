import numpy as np

from src.dirac_kahler import assemble_q6_dirac_kahler


def test_q6_dirac_kahler_is_hermitian_and_degree_odd():
    operator,audit=assemble_q6_dirac_kahler()
    assert operator.shape==(496,496)
    assert audit["cochain_dimensions"]==[64,192,240]
    assert audit["d_squared_residual"]==0.0
    assert audit["hermiticity_residual"]==0.0
    assert audit["degree_grading_anticommutator_residual"]==0.0


def test_q6_dirac_kahler_square_is_discrete_hodge_laplacian():
    _,audit=assemble_q6_dirac_kahler()
    assert audit["square_to_hodge_laplacian_residual"]<1e-12
