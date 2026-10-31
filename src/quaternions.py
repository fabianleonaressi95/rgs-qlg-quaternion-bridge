"""Quaternion algebra for the RGS-QLG Quaternion Bridge.

Conventions
-----------
- Components are stored in scalar-first order: (w, x, y, z).
- Multiplication follows the Hamilton convention.
- Quaternions represent rotations using active, right-handed rotations
  of column vectors in the companion rotations module.
- The conjugate of q is q* = (w, -x, -y, -z).
- For a unit quaternion, q^{-1} = q*.

This module implements algebraic operations only. It does not assign
a physical interpretation to quaternion components.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Iterator


DEFAULT_TOLERANCE = 1e-12


@dataclass(frozen=True, slots=True)
class Quaternion:
    """Immutable scalar-first Hamilton quaternion."""

    w: float
    x: float
    y: float
    z: float

    def __post_init__(self) -> None:
        """Reject non-finite components."""
        values = (self.w, self.x, self.y, self.z)
        if not all(math.isfinite(float(value)) for value in values):
            raise ValueError("Quaternion components must be finite.")

        # Normalize numeric inputs to Python floats.
        object.__setattr__(self, "w", float(self.w))
        object.__setattr__(self, "x", float(self.x))
        object.__setattr__(self, "y", float(self.y))
        object.__setattr__(self, "z", float(self.z))

    @classmethod
    def from_iterable(cls, values: Iterable[float]) -> Quaternion:
        """Create a quaternion from exactly four scalar-first values."""
        components = tuple(values)
        if len(components) != 4:
            raise ValueError(
                "A quaternion requires exactly four components (w, x, y, z)."
            )
        return cls(*components)

    @classmethod
    def identity(cls) -> Quaternion:
        """Return the multiplicative identity (1, 0, 0, 0)."""
        return cls(1.0, 0.0, 0.0, 0.0)

    @classmethod
    def zero(cls) -> Quaternion:
        """Return the additive identity (0, 0, 0, 0)."""
        return cls(0.0, 0.0, 0.0, 0.0)

    @classmethod
    def pure(cls, x: float, y: float, z: float) -> Quaternion:
        """Create a pure quaternion (0, x, y, z)."""
        return cls(0.0, x, y, z)

    def __iter__(self) -> Iterator[float]:
        return iter((self.w, self.x, self.y, self.z))

    def as_tuple(self) -> tuple[float, float, float, float]:
        """Return components in scalar-first order."""
        return (self.w, self.x, self.y, self.z)

    @property
    def vector(self) -> tuple[float, float, float]:
        """Return the vector part (x, y, z)."""
        return (self.x, self.y, self.z)

    @property
    def scalar(self) -> float:
        """Return the scalar part."""
        return self.w

    def __add__(self, other: Quaternion) -> Quaternion:
        if not isinstance(other, Quaternion):
            return NotImplemented
        return Quaternion(
            self.w + other.w,
            self.x + other.x,
            self.y + other.y,
            self.z + other.z,
        )

    def __sub__(self, other: Quaternion) -> Quaternion:
        if not isinstance(other, Quaternion):
            return NotImplemented
        return Quaternion(
            self.w - other.w,
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
        )

    def __neg__(self) -> Quaternion:
        return Quaternion(-self.w, -self.x, -self.y, -self.z)

    def __mul__(self, other: Quaternion | float) -> Quaternion:
        """Hamilton product, or scalar multiplication."""
        if isinstance(other, (int, float)):
            scalar = float(other)
            return Quaternion(
                scalar * self.w,
                scalar * self.x,
                scalar * self.y,
                scalar * self.z,
            )

        if not isinstance(other, Quaternion):
            return NotImplemented

        a, b, c, d = self.as_tuple()
        e, f, g, h = other.as_tuple()

        return Quaternion(
            a * e - b * f - c * g - d * h,
            a * f + b * e + c * h - d * g,
            a * g - b * h + c * e + d * f,
            a * h + b * g - c * f + d * e,
        )

    def __rmul__(self, other: float) -> Quaternion:
        if isinstance(other, (int, float)):
            return self * other
        return NotImplemented

    def conjugate(self) -> Quaternion:
        """Return q* = (w, -x, -y, -z)."""
        return Quaternion(self.w, -self.x, -self.y, -self.z)

    def norm_squared(self) -> float:
        """Return |q|^2."""
        return (
            self.w * self.w
            + self.x * self.x
            + self.y * self.y
            + self.z * self.z
        )

    def norm(self) -> float:
        """Return the Euclidean quaternion norm."""
        return math.sqrt(self.norm_squared())

    def is_unit(self, tolerance: float = DEFAULT_TOLERANCE) -> bool:
        """Check whether |q| is approximately one."""
        if tolerance < 0:
            raise ValueError("Tolerance must be non-negative.")
        return abs(self.norm() - 1.0) <= tolerance

    def normalized(self) -> Quaternion:
        """Return q / |q|. Raise for the zero quaternion."""
        magnitude = self.norm()
        if magnitude == 0.0:
            raise ZeroDivisionError("The zero quaternion cannot be normalized.")
        return Quaternion(
            self.w / magnitude,
            self.x / magnitude,
            self.y / magnitude,
            self.z / magnitude,
        )

    def inverse(self) -> Quaternion:
        """Return q^{-1} = q* / |q|^2."""
        magnitude_squared = self.norm_squared()
        if magnitude_squared == 0.0:
            raise ZeroDivisionError("The zero quaternion has no inverse.")
        conjugate = self.conjugate()
        return Quaternion(
            conjugate.w / magnitude_squared,
            conjugate.x / magnitude_squared,
            conjugate.y / magnitude_squared,
            conjugate.z / magnitude_squared,
        )

    def almost_equal(
        self,
        other: Quaternion,
        tolerance: float = DEFAULT_TOLERANCE,
    ) -> bool:
        """Componentwise approximate equality."""
        if tolerance < 0:
            raise ValueError("Tolerance must be non-negative.")
        if not isinstance(other, Quaternion):
            return False
        return all(
            abs(a - b) <= tolerance
            for a, b in zip(self.as_tuple(), other.as_tuple())
        )

    def rotate_pure(self, vector: Iterable[float]) -> Quaternion:
        """Rotate a 3-vector using q p q^{-1}.

        The quaternion q is normalized internally. The returned value
        is a pure quaternion whose vector part is the rotated vector.
        """
        values = tuple(vector)
        if len(values) != 3:
            raise ValueError("A 3D vector must contain exactly three values.")

        q = self.normalized()
        p = Quaternion.pure(*values)
        result = q * p * q.conjugate()

        # Remove only floating-point scalar leakage.
        return Quaternion(0.0, result.x, result.y, result.z)


def quaternion_dot(a: Quaternion, b: Quaternion) -> float:
    """Euclidean inner product of two quaternions."""
    return sum(x * y for x, y in zip(a.as_tuple(), b.as_tuple()))


def same_rotation(
    q1: Quaternion,
    q2: Quaternion,
    tolerance: float = DEFAULT_TOLERANCE,
) -> bool:
    """Check equality as SO(3) rotations, allowing q2 = +/- q1.

    Both inputs are normalized before comparison.
    """
    if tolerance < 0:
        raise ValueError("Tolerance must be non-negative.")

    a = q1.normalized()
    b = q2.normalized()

    direct = a.almost_equal(b, tolerance)
    antipodal = a.almost_equal(-b, tolerance)
    return direct or antipodal


def hamilton_product(
    q1: Quaternion,
    q2: Quaternion,
) -> Quaternion:
    """Explicit functional interface for the Hamilton product."""
    return q1 * q2


def quaternion_inverse(q: Quaternion) -> Quaternion:
    """Explicit functional interface for quaternion inversion."""
    return q.inverse()


def quaternion_conjugate(q: Quaternion) -> Quaternion:
    """Explicit functional interface for quaternion conjugation."""
    return q.conjugate()
