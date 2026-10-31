"""SU(2) lattice links, path transport and gauge-covariant holonomies.

Convention:
    - Links are oriented matrices U_{source,target} in SU(2).
    - Reverse links are U_{target,source} = U_{source,target}^dagger.
    - Gauge transformations act as:
          U_{gh} -> G_g U_{gh} G_h^dagger
    - Vectors are column vectors. Therefore, for a path
          v0 -> v1 -> ... -> vn
      the path transporter is
          P = U_{v(n-1),vn} ... U_{v1,v2} U_{v0,v1}.

For a closed path based at v0:
    P -> G_v0 P G_v0^dagger
and consequently Tr(P) is gauge invariant.
"""

from dataclasses import dataclass
from typing import Dict, Hashable, Iterable, List, Sequence, Tuple

import numpy as np


def dagger(matrix: np.ndarray) -> np.ndarray:
    """Hermitian conjugate."""
    return np.asarray(matrix, dtype=complex).conj().T


def is_su2(
    matrix: np.ndarray,
    atol: float = 1e-10,
) -> bool:
    """Return True if matrix satisfies the defining SU(2) conditions."""
    matrix = np.asarray(matrix, dtype=complex)

    if matrix.shape != (2, 2):
        return False

    identity = np.eye(2, dtype=complex)

    return (
        np.allclose(dagger(matrix) @ matrix, identity, atol=atol)
        and np.isclose(np.linalg.det(matrix), 1.0, atol=atol)
    )


def validate_su2(
    matrix: np.ndarray,
    atol: float = 1e-10,
) -> np.ndarray:
    """Validate and return a complex SU(2) matrix."""
    matrix = np.asarray(matrix, dtype=complex)

    if not is_su2(matrix, atol=atol):
        raise ValueError("Matrix is not an element of SU(2).")

    return matrix.copy()


@dataclass(frozen=True)
class SU2Link:
    """An oriented SU(2) link from source to target."""

    source: Hashable
    target: Hashable
    matrix: np.ndarray

    def __post_init__(self):
        object.__setattr__(
            self,
            "matrix",
            validate_su2(self.matrix),
        )

        if self.source == self.target:
            raise ValueError("A link must connect distinct vertices.")

    def reversed(self) -> "SU2Link":
        """Return the oppositely oriented link."""
        return SU2Link(
            source=self.target,
            target=self.source,
            matrix=dagger(self.matrix),
        )


class SU2Lattice:
    """Finite graph carrying oriented SU(2) link variables."""

    def __init__(self, vertices: Iterable[Hashable]):
        self.vertices = tuple(vertices)

        if len(set(self.vertices)) != len(self.vertices):
            raise ValueError("Vertex labels must be unique.")

        self._vertex_set = set(self.vertices)
        self._links: Dict[Tuple[Hashable, Hashable], np.ndarray] = {}

    def add_link(
        self,
        source: Hashable,
        target: Hashable,
        matrix: np.ndarray,
    ) -> None:
        """Add one oriented link and its inverse orientation."""
        if source not in self._vertex_set or target not in self._vertex_set:
            raise KeyError("Both link endpoints must belong to the lattice.")

        if source == target:
            raise ValueError("A link must connect distinct vertices.")

        key = (source, target)
        reverse_key = (target, source)

        if key in self._links or reverse_key in self._links:
            raise ValueError("This undirected edge already exists.")

        matrix = validate_su2(matrix)

        self._links[key] = matrix
        self._links[reverse_key] = dagger(matrix)

    def link(
        self,
        source: Hashable,
        target: Hashable,
    ) -> np.ndarray:
        """Return the oriented link matrix."""
        try:
            return self._links[(source, target)].copy()
        except KeyError as exc:
            raise KeyError(
                f"No link exists between {source!r} and {target!r} "
                "in the requested orientation."
            ) from exc

    def path_transport(
        self,
        path: Sequence[Hashable],
    ) -> np.ndarray:
        """Compute the transporter along an ordered path.

        For column vectors, later links multiply on the left:
            P = U_{v[n-1],v[n]} ... U_{v[0],v[1]}.
        """
        if len(path) < 2:
            raise ValueError("A path must contain at least two vertices.")

        if any(vertex not in self._vertex_set for vertex in path):
            raise KeyError("Every path vertex must belong to the lattice.")

        transport = np.eye(2, dtype=complex)

        for source, target in zip(path[:-1], path[1:]):
            transport = self.link(source, target) @ transport

        return transport

    def holonomy(
        self,
        closed_path: Sequence[Hashable],
    ) -> np.ndarray:
        """Compute holonomy around a closed path."""
        if len(closed_path) < 3:
            raise ValueError("A closed path must contain at least three steps.")

        if closed_path[0] != closed_path[-1]:
            raise ValueError("The path must begin and end at the same vertex.")

        return self.path_transport(closed_path)

    def wilson_trace(
        self,
        closed_path: Sequence[Hashable],
        normalized: bool = True,
    ) -> float:
        """Return the real Wilson trace of a closed path."""
        value = np.trace(self.holonomy(closed_path))

        if not np.isclose(value.imag, 0.0, atol=1e-9):
            raise ArithmeticError("SU(2) Wilson trace should be real.")

        result = float(value.real)

        return result / 2.0 if normalized else result

    def gauge_transform(
        self,
        gauge_matrices: Dict[Hashable, np.ndarray],
    ) -> "SU2Lattice":
        """Return the gauge-transformed lattice.

        Each vertex must have one SU(2) gauge matrix.
        """
        if set(gauge_matrices) != self._vertex_set:
            raise ValueError(
                "A gauge matrix must be supplied for every lattice vertex."
            )

        gauge = {
            vertex: validate_su2(gauge_matrices[vertex])
            for vertex in self.vertices
        }

        transformed = SU2Lattice(self.vertices)

        for (source, target), matrix in self._links.items():
            if self.vertices.index(source) < self.vertices.index(target):
                new_matrix = (
                    gauge[source]
                    @ matrix
                    @ dagger(gauge[target])
                )
                transformed.add_link(source, target, new_matrix)

        return transformed

    def edges(self) -> List[Tuple[Hashable, Hashable]]:
        """Return each undirected edge exactly once."""
        return [
            (source, target)
            for source, target in self._links
            if self.vertices.index(source) < self.vertices.index(target)
        ]
```
