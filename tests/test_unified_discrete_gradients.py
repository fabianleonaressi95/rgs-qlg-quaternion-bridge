"""Independent centred finite-difference checks for Q6 action derivatives."""
import numpy as np
import pytest

from src.holonomies import SU2Lattice
from src.q6_complex import build_q6_complex
from src.unified_discrete_action import ActionCouplings, UnifiedQ6Fields, unified_q6_action
from src.unified_discrete_gradients import directional_action_derivative


def _su2(axis, angle):
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    c, s = np.cos(angle / 2), np.sin(angle / 2)
    return np.array([[c-1j*s*z, -s*y-1j*s*x],
                     [s*y-1j*s*x, c+1j*s*z]], dtype=complex)


def _exp_su2(x, t):
    # su(2) is anti-Hermitian and traceless; -iX is Hermitian.
    a = float(np.sqrt(max(0.0, -0.5*np.trace(x @ x).real)))
    if a < 1e-14:
        return np.eye(2, dtype=complex) + t*x
    return np.cos(t*a)*np.eye(2) + np.sin(t*a)/a*x


def _state(seed=31):
    rng = np.random.default_rng(seed)
    q6 = build_q6_complex()
    e = np.repeat(np.eye(4)[None], 64, axis=0) + 0.025*rng.normal(size=(64,4,4))
    p = rng.normal(size=(64,2)) + 1j*rng.normal(size=(64,2))
    m = rng.normal(size=(64,2)) + 1j*rng.normal(size=(64,2))
    fields = UnifiedQ6Fields(e,p,m)
    links = SU2Lattice(q6.vertices)
    for u,v in q6.edges:
        links.add_link(u,v,_su2(rng.normal(size=3),float(rng.uniform(-0.5,0.5))))
    return rng,q6,fields,links


def _shifted(fields, links, eps, de, dp, dm, generators, q6):
    new_fields = UnifiedQ6Fields(fields.coframe+eps*de,
                                 fields.psi_plus+eps*dp,
                                 fields.psi_minus+eps*dm)
    new_links = SU2Lattice(q6.vertices)
    for u,v in q6.edges:
        U = links.link(u,v)
        X = generators.get((u,v))
        if X is not None:
            U = _exp_su2(X,eps) @ U
        new_links.add_link(u,v,U)
    return new_fields,new_links


def _fd(fields, links, q6, eps, de, dp, dm, generators, couplings):
    fp,lp = _shifted(fields,links,eps,de,dp,dm,generators,q6)
    fm,lm = _shifted(fields,links,-eps,de,dp,dm,generators,q6)
    return (unified_q6_action(fp,lp,couplings,q6)-unified_q6_action(fm,lm,couplings,q6))/(2*eps)


@pytest.mark.parametrize("sector", ["coframe","matter_plus","matter_minus","link","combined"])
def test_directional_derivatives_against_independent_central_difference(sector):
    rng,q6,fields,links = _state()
    de=np.zeros_like(fields.coframe); dp=np.zeros_like(fields.psi_plus); dm=np.zeros_like(fields.psi_minus)
    generators={}
    edge=q6.edges[37]
    X=0.17*np.array([[1j,0.4],[-0.4,-1j]],dtype=complex)
    if sector in ("coframe","combined"):
        de=0.08*rng.normal(size=de.shape)
    if sector in ("matter_plus","combined"):
        dp=0.07*(rng.normal(size=dp.shape)+1j*rng.normal(size=dp.shape))
    if sector in ("matter_minus","combined"):
        dm=0.07*(rng.normal(size=dm.shape)+1j*rng.normal(size=dm.shape))
    if sector in ("link","combined"):
        generators={edge:X}
    couplings=ActionCouplings(0.7,0.9,0.8,0.6,0.2)
    analytic=directional_action_derivative(fields,links,de,dp,dm,generators,couplings,q6)["total"]
    errors=[]
    for eps in (1e-3,3e-4,1e-4,3e-5,1e-5):
        numeric=_fd(fields,links,q6,eps,de,dp,dm,generators,couplings)
        errors.append(abs(numeric-analytic))
    assert min(errors) < 2e-6, (sector,analytic,errors)
    assert errors[-1] < 2e-5, (sector,analytic,errors)
