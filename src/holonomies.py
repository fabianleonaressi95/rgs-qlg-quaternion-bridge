"""Discrete SU(2) connections and lattice holonomies.

The module implements:
- oriented SU(2) links;
- reverse-link construction;
- path parallel transport;
- plaquette and closed-loop holonomies;
- vertex-based gauge transformations;
- gauge-invariant Wilson traces.

Convention:
    U_reverse = U_forward^{-1} = U_forward.conj().T

For a path (v0, v1, ..., vn):
    P = U[v0,v1] @ U[v1,v2] @ ... @ U[v(n-1),vn]

Gauge transformation:
    U[g,h] -> G[g] @ U[g,h] @ G[h]^{-1}
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable, Iterable, Sequence

import numpy as np
from numpy.typing import NDArray

try:
    from .su2_lifts import is_su2, su2_inverse
except ImportError:
    from su2_lifts import is_su2, su2_inverse


ComplexArray = NDArray[np.complex128]
Vertex = Hashable
Edge = tuple[Vertex, Vertex]


def _matrix2(value: NDArray[np.complex128]) -> ComplexArray:
    """Validate and copy a finite 2x2 complex matrix."""
    matrix = np.asarray(value, dtype=np.complex128)

    if matrix.shape != (2, 2):
        raise ValueError("An SU(2) link must be a 2x2 matrix.")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Matrix entries must be finite.")

    return matrix.copy()


@dataclass(frozen=True)
class SU2Link:
    """An oriented link carrying an SU(2) matrix."""

    source: Vertex
    target: Vertex
    matrix: ComplexArray

    def __post_init__(self) -> None:
        matrix = _matrix2(self.matrix)

        if self.source == self.target:
            raise ValueError("A link must connect two distinct vertices.")

        if not is_su2(matrix):
            raise ValueError("Link matrix is not in SU(2).")

        matrix.setflags(write=False)
        object.__setattr__(self, "matrix", matrix)

    def reverse(self) -> SU2Link:
        """Return the oppositely oriented link."""
        return SU2Link(
            source=self.target,
            target=self.source,
            matrix=su2_inverse(self.matrix),
        )

    def gauge_transform(
        self,
        source_gauge: ComplexArray,
        target_gauge: ComplexArray,
    ) -> SU2Link:
        """Apply U' = G_source U G_target^{-1}."""
        gs = _matrix2(source_gauge)
        gt = _matrix2(target_gauge)

        if not is_su2(gs) or not is_su2(gt):
            raise ValueError("Gauge transformations must belong to SU(2).")

        transformed = gs @ self.matrix @ gt.conj().T

        return SU2Link(
            source=self.source,
            target=self.target,
            matrix=transformed,
        )


class SU2Lattice:
    """Finite graph equipped with oriented SU(2) link variables.

    Each undirected edge is stored once in a canonical orientation.
    Reverse links are generated using the group inverse.
    """

    def __init__(self) -> None:
        self._links: dict[Edge, SU2Link] = {}
        self._vertices: set[Vertex] = set()

    @property
    def vertices(self) -> frozenset[Vertex]:
        return frozenset(self._vertices)

    @property
    def edges(self) -> tuple[Edge, ...]:
        return tuple(self._links.keys())

    def add_link(
        self,
        source: Vertex,
        target: Vertex,
        matrix: ComplexArray,
    ) -> None:
        """Add one oriented link, rejecting duplicate undirected edges."""
        if source == target:
            raise ValueError("Self-links are not supported.")

        key = (source, target)
        reverse_key = (target, source)

        if key in self._links or reverse_key in self._links:
            raise ValueError("This undirected edge is already present.")

        link = SU2Link(source, target, matrix)
        self._links[key] = link
        self._vertices.update((source, target))

    def link(self, source: Vertex, target: Vertex) -> ComplexArray:
        """Return the matrix for either orientation of an existing edge."""
        key = (source, target)

        if key in self._links:
            return self._links[key].matrix.copy()

        reverse_key = (target, source)

        if reverse_key in self._links:
            return self._links[reverse_key].reverse().matrix.copy()

        raise KeyError(f"No link connects {source!r} and {target!r}.")

    def path_transport(
        self,
        path: Sequence[Vertex],
    ) -> ComplexArray:
        """Compute ordered parallel transport along a vertex path.

        The rightmost matrix acts first on a state vector.
        """
        if len(path) < 2:
            raise ValueError("A path must contain at least two vertices.")

        transport = np.eye(2, dtype=np.complex128)

        for source, target in zip(path[:-1], path[1:]):
            transport = transport @ self.link(source, target)

        return transport

    def holonomy(
        self,
        closed_path: Sequence[Vertex],
    ) -> ComplexArray:
        """Compute the holonomy of a closed path."""
        if len(closed_path) < 3:
            raise ValueError("A closed path must contain at least three entries.")

        if closed_path[0] != closed_path[-1]:
            raise ValueError("The path is not closed.")

        return self.path_transport(closed_path)

    def plaquette_holonomy(
        self,
        vertices: Sequence[Vertex],
    ) -> ComplexArray:
        """Compute a four-edge plaquette holonomy.

        Expected input:
            (v0, v1, v2, v3, v0)
        """
        if len(vertices) != 5:
            raise ValueError(
                "A plaquette requires four vertices and a repeated endpoint."
            )

        return self.holonomy(vertices)

    def wilson_trace(
        self,
        closed_path: Sequence[Vertex],
        normalized: bool = True,
    ) -> complex:
        """Return the trace of a closed-loop holonomy.

        For the fundamental SU(2) representation, normalized Wilson trace
        is Tr(H)/2.
        """
        h = self.holonomy(closed_path)
        value = complex(np.trace(h))

        return value / 2.0 if normalized else value

    def gauge_transform(
        self,
        gauge: dict[Vertex, ComplexArray],
    ) -> SU2Lattice:
        """Return a new lattice transformed by vertex SU(2) matrices."""
        missing = self._vertices.difference(gauge.keys())
        if missing:
            raise ValueError(f"Missing gauge matrices for vertices: {missing!r}")

        for vertex in self._vertices:
            if not is_su2(gauge[vertex]):
                raise ValueError(
                    f"Gauge matrix at vertex {vertex!r} is not in SU(2)."
                )

        transformed = SU2Lattice()

        for (source, target), link in self._links.items():
            transformed_link = link.gauge_transform(
                gauge[source],
                gauge[target],
            )
            transformed.add_link(
                source,
                target,
                transformed_link.matrix,
            )

        return transformed

    def reverse_link_residual(
        self,
        source: Vertex,
        target: Vertex,
    ) -> float:
        """Check U[target,source] = U[source,target] dagger."""
        forward = self.link(source, target)
        reverse = self.link(target, source)

        return float(
            np.linalg.norm(reverse - forward.conj().T, ord="fro")
        )


def wilson_trace(
    holonomy: ComplexArray,
    normalized: bool = True,
) -> complex:
    """Compute the trace of an already constructed SU(2) holonomy."""
    h = _matrix2(holonomy)

    if not is_su2(h):
        raise ValueError("Holonomy matrix is not in SU(2).")

    value = complex(np.trace(h))
    return value / 2.0 if normalized else value


def gauge_covariance_residual(
    original: ComplexArray,
    transformed: ComplexArray,
    base_gauge: ComplexArray,
) -> float:
    """Check H' = G H G^{-1} for a based closed-loop holonomy."""
    h = _matrix2(original)
    hp = _matrix2(transformed)
    g = _matrix2(base_gauge)

    if not is_su2(h) or not is_su2(hp) or not is_su2(g):
        raise ValueError("All matrices must belong to SU(2).")

    expected = g @ h @ g.conj().T

    return float(np.linalg.norm(hp - expected, ord="fro"))


def wilson_gauge_invariance_residual(
    original: ComplexArray,
    transformed: ComplexArray,
) -> float:
    """Check invariance of the normalized Wilson trace."""
    w1 = wilson_trace(original)
    w2 = wilson_trace(transformed)

    return float(abs(w1 - w2))


def holonomy_unitarity_residual(
    holonomy: ComplexArray,
) -> float:
    """Return ||H†H-I||_F."""
    h = _matrix2(holonomy)
    return float(
        np.linalg.norm(
            h.conj().T @ h - np.eye(2),
            ord="fro",
        )
    )


def holonomy_determinant_residual(
    holonomy: ComplexArray,
) -> float:
    """Return |det(H)-1|."""
    h = _matrix2(holonomy)
    return float(abs(np.linalg.det(h) - 1.0))
