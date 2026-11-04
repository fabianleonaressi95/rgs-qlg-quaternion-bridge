import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist

from src.rgs_psu_cut_project import build_cut_project_data


def test_rhombic_triacontahedral_window_has_30_rhombic_facets():
    data = build_cut_project_data()
    hull = ConvexHull(data.q6_internal)
    planes = np.unique(np.round(hull.equations, 9), axis=0)
    assert len(planes) == 30
    edge_spreads = []
    for plane in planes:
        ids = np.flatnonzero(
            np.abs(data.window_vertices @ plane[:3] + plane[3]) < 1e-8
        )
        assert len(ids) == 4
        edge_lengths = np.sort(pdist(data.window_vertices[ids]))[:4]
        edge_spreads.append(np.ptp(edge_lengths))
    assert max(edge_spreads) < 1e-10


def test_q6_edges_project_to_equal_length_icosahedral_vectors():
    data = build_cut_project_data()
    edges = []
    for vertex in range(64):
        for axis in range(6):
            other = vertex ^ (1 << (5 - axis))
            if vertex < other:
                edges.append((vertex, other))
    edges = np.asarray(edges, dtype=np.int64)
    lengths = np.linalg.norm(
        data.q6_parallel[edges[:, 1]] - data.q6_parallel[edges[:, 0]], axis=1
    )
    assert edges.shape == (192, 2)
    assert np.allclose(lengths, 1.0 / np.sqrt(2.0), atol=1e-14)
