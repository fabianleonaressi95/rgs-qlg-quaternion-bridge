"""Golden ratio geometry and supergolden constant.

This module implements:
- Golden ratio φ = (1 + √5)/2
- Supergolden constant ψ (where ψ³ = ψ + 1)
- Golden icosahedron coordinates
- Star tetrahedron geometry
"""

from .constants import PHI, PSI, PHI_INVERSE, PHI_SQUARED
from .icosahedron import icosahedron_vertices, dodecahedron_vertices
from .star_tetrahedron import star_tetrahedron_vertices

__all__ = [
    "PHI",
    "PSI",
    "PHI_INVERSE",
    "PHI_SQUARED",
    "icosahedron_vertices",
    "dodecahedron_vertices",
    "star_tetrahedron_vertices",
]
