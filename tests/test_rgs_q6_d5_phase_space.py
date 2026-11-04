import numpy as np
import pytest
from src.rgs_q6_d5_phase_space import (
 build_q6_d5_phase_space_sector,d5_type_iv_contains,phase_space_to_d5_tangent
)


def test_q6_top_coexact_sector_is_five_dimensional_A5_irrep():
 r=build_q6_d5_phase_space_sector()
 assert r["quotient_configuration_dimension"]==129
 assert r["selected_mode_multiplicity"]==5
 assert r["A5_irrep_multiplicities"]=={"1":0,"3":0,"3_prime":0,"4":0,"5":1}
 assert r["all_checks_passed"]


def test_cotangent_lift_is_ten_dimensional_kahler_and_symplectic():
 r=build_q6_d5_phase_space_sector(beta=2.4)
 assert r["phase_space_dimension"]==10
 assert r["checks"]["complex_structure"]
 assert r["checks"]["symplectic_compatibility"]
 assert r["checks"]["kahler_metric_positive"]
 assert r["checks"]["SO2_phase_isometry"]
 assert r["checks"]["A5_phase_isometry"]


def test_type_iv_domain_contains_origin_and_small_tangent_vectors():
 assert d5_type_iv_contains(np.zeros(5,dtype=complex))
 z=phase_space_to_d5_tangent(np.full(5,1e-3),np.full(5,-2e-3))
 assert d5_type_iv_contains(z)
 assert not d5_type_iv_contains(np.full(5,2.,dtype=complex))


def test_phase_space_to_d5_tangent_has_expected_complex_dimension():
 q=np.array([.01,.02,0.,-.01,.005])
 p=np.array([0.,.01,-.02,.005,.01])
 z=phase_space_to_d5_tangent(q,p,beta=1.3)
 assert z.shape==(5,)
 assert np.iscomplexobj(z)
 assert np.all(np.isfinite(z))


@pytest.mark.parametrize("bad", [0.,-1.,np.inf,np.nan])
def test_invalid_beta_is_rejected(bad):
 with pytest.raises(ValueError):
  build_q6_d5_phase_space_sector(beta=bad)
