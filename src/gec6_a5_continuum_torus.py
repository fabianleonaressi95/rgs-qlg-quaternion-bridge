"""Exact spectral data for the flat A5 root-lattice torus."""
from itertools import product
from math import pi

Q = (
    (5, 4, 3, 2, 1),
    (4, 8, 6, 4, 2),
    (3, 6, 9, 6, 3),
    (2, 4, 6, 8, 4),
    (1, 2, 3, 4, 5),
)


def quadratic_form(k):
    if len(k) != 5 or any(not isinstance(x, int) for x in k):
        raise ValueError("k must be a five-vector of integers")
    return sum(k[i] * Q[i][j] * k[j] for i in range(5) for j in range(5))


def shortest_dual_vectors():
    # Q > (3/2)I implies q(k)<=5 => ||k||^2<10/3.
    # Thus {-1,0,1}^5 is an exhaustive finite search set.
    return tuple(k for k in product((-1, 0, 1), repeat=5)
                 if k != (0, 0, 0, 0, 0) and quadratic_form(k) == 5)


def shortest_vector_certificate():
    values = [quadratic_form(k) for k in product((-1, 0, 1), repeat=5)
              if k != (0, 0, 0, 0, 0)]
    vectors = shortest_dual_vectors()
    return {
        "leading_principal_minors_Q_minus_3_over_2_I":
            ["7/2", "27/4", "81/8", "189/16", "243/32"],
        "finite_search_set": "{-1,0,1}^5",
        "candidate_count_including_zero": 243,
        "minimum": min(values),
        "minimizers": [list(k) for k in vectors],
        "multiplicity": len(vectors),
    }


def continuum_eigenvalue(k):
    return (2.0 * pi**2 / 3.0) * quadratic_form(k)


def continuum_spectral_data():
    cert = shortest_vector_certificate()
    if cert["minimum"] != 5 or cert["multiplicity"] != 12:
        raise AssertionError("shortest dual-vector certificate failed")
    return {
        "space": "T_A5=V/A5, V={x in R^6: sum_i x_i=0}",
        "dimension": 5,
        "dual_lattice": "A5*",
        "eigenvalue_formula": "lambda(k)=(2*pi^2/3) k^T Q k",
        "first_quadratic_value": 5,
        "first_eigenvalue": "10*pi^2/3",
        "first_eigenvalue_numeric": 10.0*pi**2/3.0,
        "first_eigenvalue_multiplicity": 12,
        "shortest_vector_certificate": cert,
    }


def de_aquino_comparison():
    return {
        "mass_ratio": "m_g/m_i0 = 3 - 2*sqrt(1 + xi^2)",
        "zero_threshold_squared": "5/4",
        "a5_shortest_quadratic_value": 5,
        "candidate_normalization": "xi_k^2 = (k^T Q k)/4",
        "formal_match": "xi_k^2=5/4 on the twelve shortest dual vectors",
        "status": "CONJECTURAL NORMALIZATION ONLY; no physical momentum map or operator intertwiner is established",
        "warning": "The earlier artificial A6/Dirac grid audit is a separate numerical experiment, not a proof of this bridge.",
    }


if __name__ == "__main__":
    import json
    data = continuum_spectral_data()
    data["de_aquino_comparison"] = de_aquino_comparison()
    print(json.dumps(data, indent=2, default=str))
