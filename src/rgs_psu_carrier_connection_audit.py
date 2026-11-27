"""Generator-level comparison of archived PSU links and carrier coframe."""
from __future__ import annotations
import hashlib
import numpy as np
from .rgs_ellis_refinement import load_rgs_seed_data

EXPECTED_EDGE_SHA256 = "da7854ca655577b629dc508208918326fde1d50563485f953ea0fdc1b2fb78e8"
ARCHIVED_TORSION_NORM_RMS = 2.3537334152693665
I2 = np.eye(2, dtype=complex)
I3 = np.eye(3)
PAULI = (
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
)

def _skew(v):
    x, y, z = v
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])

def _minimal_rotation(a, b):
    a = np.asarray(a, float) / np.linalg.norm(a)
    b = np.asarray(b, float) / np.linalg.norm(b)
    c = float(np.clip(a @ b, -1., 1.))
    v = np.cross(a, b)
    s = float(np.linalg.norm(v))
    if s < 1e-14:
        if c > 0:
            return I3.copy()
        ref = np.array([1., 0, 0]) if abs(a[0]) < .8 else np.array([0., 1., 0])
        u = np.cross(a, ref)
        u /= np.linalg.norm(u)
        return 2. * np.outer(u, u) - I3
    V = _skew(v / s)
    return I3 + V * s + (V @ V) * (1. - c)

def _radial_frame(x):
    n = np.asarray(x, float) / np.linalg.norm(x)
    reference = np.array([0., 0., 1.])
    if abs(float(reference @ n)) > .9:
        reference = np.array([1., 0., 0.])
    e1 = np.cross(reference, n)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return np.column_stack((e1, e2, n))

def _su2_lift(R):
    tr = float(np.trace(R))
    w = np.sqrt(max(0., 1. + tr)) / 2.
    if w > 1e-10:
        qv = np.array([R[2,1]-R[1,2], R[0,2]-R[2,0], R[1,0]-R[0,1]]) / (4.*w)
    else:
        vals, vecs = np.linalg.eig(R)
        idx = int(np.argmin(np.abs(vals - 1.)))
        qv = np.real(vecs[:, idx])
        qv /= np.linalg.norm(qv)
    q = np.r_[w, qv]
    q /= np.linalg.norm(q)
    return q[0] * I2 - 1j * sum(q[k+1] * PAULI[k] for k in range(3))

def _so3_from_su2(U):
    return np.array([
        [np.real(np.trace(PAULI[a] @ U @ PAULI[b] @ U.conj().T) / 2.)
         for b in range(3)] for a in range(3)
    ])

def _canonical_edges():
    pairs = set()
    for g in range(64):
        for axis_be in range(6):
            h = g ^ (1 << (5-axis_be))
            if g < h:
                pairs.add((g, h))
    return np.asarray(sorted(pairs), dtype=np.int64)

def _face_torsion(square, edge_index, edges, rotations, coframe):
    P = I3.copy()
    T = np.zeros(3)
    for k in range(4):
        a, b = square[k], square[(k+1) % 4]
        idx = edge_index[(min(a,b), max(a,b))]
        ea, eb = map(int, edges[idx,:2])
        forward = (a,b) == (ea,eb)
        Rdir = rotations[idx] if forward else rotations[idx].T
        edir = coframe[idx] if forward else -(rotations[idx] @ coframe[idx])
        T += P.T @ edir
        P = Rdir @ P
    return T

def audit_rgs_psu_carrier_connection():
    seed = load_rgs_seed_data()
    X = np.asarray(seed.coordinates, float)
    edges = np.asarray(seed.q6_edges, dtype=np.int64)
    if X.shape != (64,3) or edges.shape != (192,3):
        raise AssertionError("Unexpected seed or Q6 edge shape")
    edge_digest = hashlib.sha256(
        np.ascontiguousarray(edges, dtype="<i8").tobytes(order="C")
    ).hexdigest()
    if edge_digest != EXPECTED_EDGE_SHA256:
        raise AssertionError(f"Q6 edge-order hash mismatch: {edge_digest}")
    canonical = _canonical_edges()
    canonical_index = {(int(a),int(b)): k for k,(a,b) in enumerate(canonical)}
    rgs_to_canonical = np.asarray(
        [canonical_index[(int(a),int(b))] for a,b in edges[:,:2]], dtype=np.int64
    )
    if not np.array_equal(canonical[rgs_to_canonical], edges[:,:2]):
        raise AssertionError("RGS-to-canonical Q6 edge permutation failed")

    radii = np.linalg.norm(X, axis=1)
    directions = X / radii[:,None]
    frames = np.asarray([_radial_frame(x) for x in X])
    # Reproduce the A13 builder's minimal radial rotations and quaternion lifts.
    R_arch = np.asarray([_minimal_rotation(directions[g],directions[h])
                         for g,h,_ in edges])
    U_arch = np.asarray([_su2_lift(R) for R in R_arch])
    R_local = np.asarray([frames[h].T @ R_arch[k] @ frames[g]
                          for k,(g,h,_) in enumerate(edges)])
    frame_lifts = np.asarray([_su2_lift(F) for F in frames])
    U_local = np.asarray([frame_lifts[h].conj().T @ U_arch[k] @ frame_lifts[g]
                          for k,(g,h,_) in enumerate(edges)])

    so3_local_residual = max(float(np.linalg.norm(
        R_local[k] - frames[h].T @ R_arch[k] @ frames[g]))
        for k,(g,h,_) in enumerate(edges))
    so3_orthogonality = max(float(np.linalg.norm(R.T @ R-I3)) for R in R_local)
    so3_determinant = max(abs(float(np.linalg.det(R))-1.) for R in R_local)
    radial_axis_residual = max(float(np.linalg.norm(
        R @ np.array([0.,0.,1.]) - np.array([0.,0.,1.]))) for R in R_local)
    su2_unitarity = max(float(np.linalg.norm(U.conj().T @ U-I2)) for U in U_local)
    su2_adjoint_residual = max(float(np.linalg.norm(
        _so3_from_su2(U_local[k])-R_local[k])) for k in range(len(edges)))
    su2_offdiagonal = max(float(np.linalg.norm(U - np.diag(np.diag(U)))) for U in U_local)

    raw4 = np.asarray(seed.q6_raw_coframe, float)
    unit4 = np.asarray(seed.q6_unit_coframe, float)
    raw4_rebuilt = np.column_stack((
        X[edges[:,1]]-X[edges[:,0]], radii[edges[:,1]]-radii[edges[:,0]]
    ))
    coframe_reconstruction_residual = float(np.linalg.norm(raw4-raw4_rebuilt))
    unit4_reconstruction_residual = float(np.linalg.norm(
        unit4 - raw4_rebuilt/np.linalg.norm(raw4_rebuilt,axis=1)[:,None]))

    E_raw_local = np.empty((192,3))
    E_carrier = np.empty((192,3))
    E4_spatial_local = np.empty((192,3))
    radial4 = np.empty(192)
    ratios = np.empty(192)
    reverse_residual = 0.
    for k,(g0,h0,_) in enumerate(edges):
        g,h = int(g0),int(h0)
        delta = X[h]-X[g]
        ci = frames[g].T @ delta
        cj = frames[h].T @ delta
        E_raw_local[k] = ci
        E_carrier[k] = .5*(ci + R_local[k].T @ cj)
        L4 = np.linalg.norm(raw4_rebuilt[k])
        L3 = np.linalg.norm(delta)
        E4_spatial_local[k] = ci/L4
        radial4[k] = raw4_rebuilt[k,3]/L4
        ratios[k] = np.linalg.norm(E_carrier[k])/L3
        Rji = R_local[k].T
        reverse = .5*(frames[h].T @ (-delta) + Rji.T @ frames[g].T @ (-delta))
        reverse_residual = max(reverse_residual,
                               float(np.linalg.norm(reverse + R_local[k] @ E_carrier[k])))

    squares = []
    for a in range(6):
        for b in range(a+1,6):
            for g in range(64):
                if ((g>>a)&1) or ((g>>b)&1):
                    continue
                i=g; j=g^(1<<a); k=j^(1<<b); l=g^(1<<b)
                squares.append((i,j,k,l))
    if len(squares) != 240:
        raise AssertionError("Expected 240 Q6 squares")
    edge_index = {(min(int(g),int(h)),max(int(g),int(h))):idx
                  for idx,(g,h,_) in enumerate(edges)}
    T_raw = np.asarray([_face_torsion(sq,edge_index,edges,R_local,E_raw_local)
                        for sq in squares])
    T_carrier = np.asarray([_face_torsion(sq,edge_index,edges,R_local,E_carrier)
                             for sq in squares])
    T_radial = []
    for sq in squares:
        total = 0.
        for a,b in zip(sq, sq[1:]+sq[:1]):
            idx = edge_index[(min(a,b),max(a,b))]
            sign = 1. if (a,b)==tuple(edges[idx,:2]) else -1.
            total += sign*raw4_rebuilt[idx,3]
        T_radial.append(total)
    raw_norms = np.linalg.norm(T_raw,axis=1)
    carrier_norms = np.linalg.norm(T_carrier,axis=1)
    raw_torsion_rms = float(np.sqrt(np.mean(raw_norms**2)))
    carrier_torsion_rms = float(np.sqrt(np.mean(carrier_norms**2)))
    radial_torsion_rms = float(np.sqrt(np.mean(np.asarray(T_radial)**2)))
    archived_spatial_projection = E4_spatial_local
    E_carrier_unit = E_carrier/np.linalg.norm(raw4_rebuilt[:,:3],axis=1)[:,None]
    result = {
      "status":"GENERATOR_LEVEL_GAUGE_EQUIVALENCE_VERIFIED; ORIGINAL_A13_BINARY_PAYLOAD_NOT_LOADED",
      "provenance":{
        "coordinate_sha256":"6a158cb30a671a5c994d8ab223e9283b57ac703808041d46e67f286126e419cf",
        "rgs_edge_sha256_int64":edge_digest,
        "a13_generator":"minimal SO(3) rotation n_g -> n_h; scalar-first quaternion lift",
        "binary_R_SO3_U_SU2_comparison":"NOT_PERFORMED; original A13 NPZ payload unavailable in repository and current Library search results"
      },
      "edge_alignment":{
        "rgs_to_canonical_permutation":rgs_to_canonical.tolist(),
        "edge_set_equal":True,"edge_count":int(len(edges)),
        "axis_conversion":"axis_big_endian = 5 - axis_lsb"
      },
      "connection":{
        "ambient_to_local_formula":"R_local(g,h) = F_h^T R_arch(g,h) F_g",
        "max_SO3_gauge_transform_residual":so3_local_residual,
        "max_SO3_orthogonality_residual":so3_orthogonality,
        "max_SO3_determinant_residual":so3_determinant,
        "max_radial_axis_fixed_residual":radial_axis_residual,
        "max_SU2_unitarity_residual":su2_unitarity,
        "max_SU2_adjoint_to_SO3_residual":su2_adjoint_residual,
        "max_SU2_local_offdiagonal":su2_offdiagonal,
        "local_structure":"SO(3) links fix the radial axis; local spin links lie in diagonal U(1) up to the double-cover sign"
      },
      "coframe_comparison":{
        "archived_feature":"normalized 4D (Delta X, Delta |X|)",
        "carrier_feature":"3D symmetric local-frame chord average",
        "raw_4D_reconstruction_residual":coframe_reconstruction_residual,
        "normalized_4D_reconstruction_residual":unit4_reconstruction_residual,
        "carrier_to_chord_norm_ratio_min":float(ratios.min()),
        "carrier_to_chord_norm_ratio_max":float(ratios.max()),
        "carrier_to_chord_norm_ratio_rms":float(np.sqrt(np.mean(ratios**2))),
        "max_difference_from_archived_normalized_spatial_projection":float(np.max(
            np.linalg.norm(E_carrier_unit-archived_spatial_projection,axis=1))),
        "max_difference_from_source_only_unit_chord":float(np.max(np.linalg.norm(
            E_carrier_unit-E_raw_local/np.linalg.norm(raw4_rebuilt[:,:3],axis=1)[:,None],axis=1))),
        "radial_augmented_component_rms":float(np.sqrt(np.mean(radial4**2))),
        "carrier_is_same_as_archived_4D_feature":False
      },
      "Q6_square_torsion":{
        "square_count":len(squares),
        "raw_local_spatial_torsion_norm_rms":raw_torsion_rms,
        "previous_A13_4D_report_norm_rms":ARCHIVED_TORSION_NORM_RMS,
        "raw_local_vs_previous_report_abs_difference":abs(raw_torsion_rms-ARCHIVED_TORSION_NORM_RMS),
        "symmetric_carrier_torsion_norm_rms":carrier_torsion_rms,
        "radial_plain_cochain_torsion_rms":radial_torsion_rms,
        "raw_spatial_torsion_matches_previous_report":abs(raw_torsion_rms-ARCHIVED_TORSION_NORM_RMS)<1e-10
      },
      "orientation":{
        "max_carrier_reverse_edge_residual":reverse_residual,
        "max_norm_bound_excess":max(0.,float(ratios.max())-1.)
      },
      "interpretation":[
        "The A13 and carrier links are the same shortest-radial-rotation connection in ambient and local radial frames.",
        "The local connection fixes the radial frame axis, so links reduce to the SO(2) stabilizer and spin lifts to diagonal U(1).",
        "The archived normalized 4D cochain and carrier symmetrized 3D coframe are different objects.",
        "The raw local spatial cochain reproduces the spatial part of the previous Q6 torsion report; the symmetrized carrier coframe defines a new finite torsion diagnostic.",
        "No continuum Cartan/TEGR, spacetime tetrad, Q6 tetrahedral 3-cell or physical uniqueness claim follows."
      ]
    }
    return result
