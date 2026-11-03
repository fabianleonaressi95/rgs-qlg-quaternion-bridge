"""Complete audit of quaternion bridge and 2I proximity.

This script:
1. Loads RGS vierbein from covariance matrices
2. Constructs SU(2) connection M_ij = E_j * E_i^{-1}
3. Performs polar decomposition R_ij * S_ij
4. Lifts to quaternion link variables
5. Computes distance from binary icosahedral group 2I
6. Calculates holonomies on all 240 Q6 plaquettes
7. Analyzes curvature distribution
8. Generates comprehensive audit report
"""

import sys
import numpy as np
from pathlib import Path
import json

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent))

from vierbein.construction import build_vierbein
from quaternion_bridge import (
    su2_connection_from_vierbein,
    analyze_2I_proximity,
    compute_all_holonomies,
    analyze_holonomy_distribution,
)
from binary_icosahedral import binary_icosahedral_elements
from q6 import q6_edges, q6_plaquettes
from load_rgs_matrices import load_rgs_matrices


def main():
    """Main execution function."""
    print("=" * 70)
    print("RGS-PSU QUATERNION BRIDGE AND 2I PROXIMITY AUDIT")
    print("=" * 70)
    print()

    # Load RGS matrices
    print("Loading RGS matrices...")
    matrices, metadata = load_rgs_matrices()
    print(f"Loaded {len(matrices)} arrays")
    print(f"SHA-256: {metadata['sha256']}")
    print()

    # Extract vertex positions
    if "X_Oh_reconstructed" not in matrices:
        print("ERROR: X_Oh_reconstructed not found")
        sys.exit(1)

    X = matrices["X_Oh_reconstructed"]
    print(f"Vertex positions: shape {X.shape}")
    print()

    # Historical parameters
    ell_squared = 7.035514834154695
    c0 = 43.96495130838031

    # Build vierbein
    print("STEP 1: VIERBEIN CONSTRUCTION")
    print("-" * 70)
    E, vierbein_metadata = build_vierbein(X, ell_squared, c0)
    print(f"Vierbein shape: {E.shape}")
    print(f"Positive definite: {vierbein_metadata['positive_definite_count']}/{vierbein_metadata['total_matrices']}")
    print()

    # Get Q6 adjacency
    print("STEP 2: Q6 STRUCTURE")
    print("-" * 70)
    adjacency = q6_edges()
    print(f"Q6 edges: {len(adjacency)}")
    print(f"Expected: 192")
    print(f"Correct: {len(adjacency) == 192}")
    print()

    # Build SU(2) connection
    print("STEP 3: SU(2) CONNECTION CONSTRUCTION")
    print("-" * 70)
    print("Computing M_ij = E_j * E_i^{-1}")
    print("Polar decomposition: M_ij = R_ij * S_ij")
    print("Lifting to SU(2): R_ij → q_ij")
    print()

    su2_connection = su2_connection_from_vierbein(E, adjacency)
    print(f"SU(2) link variables: {len(su2_connection)}")
    print(f"Expected: 192")
    print(f"Correct: {len(su2_connection) == 192}")
    print()

    # Load 2I elements
    print("STEP 4: BINARY ICOSAHEDRAL GROUP 2I")
    print("-" * 70)
    elements_2I = binary_icosahedral_elements()
    print(f"2I elements: {elements_2I.shape[0]}")
    print(f"Expected: 120")
    print(f"Correct: {elements_2I.shape[0] == 120}")
    print()

    # Analyze 2I proximity
    print("STEP 5: 2I PROXIMITY ANALYSIS")
    print("-" * 70)
    proximity_analysis = analyze_2I_proximity(su2_connection, elements_2I)

    print(f"Number of links: {proximity_analysis['num_links']}")
    print(f"Mean distance from 2I: {proximity_analysis['mean_distance']:.6f}")
    print(f"Std distance: {proximity_analysis['std_distance']:.6f}")
    print(f"Min distance: {proximity_analysis['min_distance']:.6f}")
    print(f"Max distance: {proximity_analysis['max_distance']:.6f}")
    print(f"Median distance: {proximity_analysis['median_distance']:.6f}")
    print(f"Fraction close (< 0.1): {proximity_analysis['fraction_close']:.2%}")
    print()

    # Interpret proximity
    mean_dist = proximity_analysis['mean_distance']
    if mean_dist < 0.1:
        print("✓ STRONG PROXIMITY: Connection is very close to 2I")
        print("  → The geometric connection spontaneously approaches the binary icosaededral structure")
    elif mean_dist < 0.5:
        print("⚠ MODERATE PROXIMITY: Connection shows some affinity to 2I")
        print("  → Partial alignment with 2I structure")
    else:
        print("✗ WEAK PROXIMITY: Connection is not close to 2I")
        print("  → The geometric connection does not spontaneously approach 2I")
    print()

    # Compute holonomies
    print("STEP 6: PLAQUETTE HOLONOMIES")
    print("-" * 70)
    plaquettes = q6_plaquettes()
    print(f"Q6 plaquettes: {len(plaquettes)}")
    print(f"Expected: 240")
    print(f"Correct: {len(plaquettes) == 240}")
    print()

    print("Computing holonomies U_f = ∏ q_e^{ε_fe}...")
    holonomies = compute_all_holonomies(su2_connection, plaquettes)
    print(f"Holonomies computed: {len(holonomies)}")
    print()

    # Analyze holonomy distribution
    print("STEP 7: HOLONOMY DISTRIBUTION ANALYSIS")
    print("-" * 70)
    holonomy_analysis = analyze_holonomy_distribution(holonomies, elements_2I)

    ca = holonomy_analysis['curvature_angles']
    print(f"Curvature angles:")
    print(f"  Mean: {ca['mean']:.6f} rad ({np.degrees(ca['mean']):.2f}°)")
    print(f"  Std: {ca['std']:.6f} rad ({np.degrees(ca['std']):.2f}°)")
    print(f"  Min: {ca['min']:.6f} rad ({np.degrees(ca['min']):.2f}°)")
    print(f"  Max: {ca['max']:.6f} rad ({np.degrees(ca['max']):.2f}°)")
    print(f"  Median: {ca['median']:.6f} rad ({np.degrees(ca['median']):.2f}°)")
    print()

    print(f"Flat plaquettes (angle ≈ 0): {holonomy_analysis['flat_plaquettes']}/{holonomy_analysis['num_plaquettes']}")
    print()

    d2I = holonomy_analysis['2I_distances']
    print(f"Holonomy distances from 2I:")
    print(f"  Mean: {d2I['mean']:.6f}")
    print(f"  Std: {d2I['std']:.6f}")
    print(f"  Min: {d2I['min']:.6f}")
    print(f"  Max: {d2I['max']:.6f}")
    print(f"  Median: {d2I['median']:.6f}")
    print()

    # Interpret holonomy results
    if holonomy_analysis['flat_plaquettes'] > 0:
        print("⚠ Some plaquettes are flat (zero curvature)")
        print("  → May indicate degenerate geometry or special symmetries")
    else:
        print("✓ All plaquettes have non-zero curvature")
        print("  → Genuine geometric curvature present")
    print()

    # Comparison with link-level 2I proximity
    link_mean = proximity_analysis['mean_distance']
    holonomy_mean = d2I['mean']

    print("STEP 8: LINK vs HOLONOMY 2I PROXIMITY")
    print("-" * 70)
    print(f"Link-level mean distance: {link_mean:.6f}")
    print(f"Holonomy-level mean distance: {holonomy_mean:.6f}")
    print(f"Ratio: {holonomy_mean / link_mean:.3f}")
    print()

    if holonomy_mean < link_mean:
        print("✓ Holonomies are CLOSER to 2I than individual links")
        print("  → Curvature operations enhance 2I alignment")
    elif holonomy_mean > link_mean:
        print("✗ Holonomies are FARTHER from 2I than individual links")
        print("  → Curvature operations reduce 2I alignment")
    else:
        print("≈ Holonomies have SIMILAR 2I proximity to links")
        print("  → Curvature preserves 2I alignment")
    print()

    # Save results
    print("STEP 9: SAVING RESULTS")
    print("-" * 70)

    results_dir = Path("results/json")
    results_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "metadata": {
            "file_hash": metadata["sha256"],
            "ell_squared": ell_squared,
            "c0": c0,
            "num_vertices": int(X.shape[0]),
        },
        "vierbein": vierbein_metadata,
        "connection": {
            "num_links": len(su2_connection),
            "num_edges_expected": 192,
        },
        "2I_proximity": proximity_analysis,
        "holonomies": {
            "num_plaquettes": len(holonomies),
            "num_expected": 240,
        },
        "holonomy_analysis": holonomy_analysis,
        "interpretation": {
            "link_2I_strong": mean_dist < 0.1,
            "link_2I_moderate": mean_dist < 0.5,
            "holonomy_enhances_2I": holonomy_mean < link_mean,
            "flat_plaquettes_exist": holonomy_analysis['flat_plaquettes'] > 0,
        },
    }

    # Convert numpy types
    def convert_types(obj):
        if isinstance(obj, (np.integer, int)):
            return int(obj)
        elif isinstance(obj, (np.floating, float)):
            return float(obj)
        elif isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_types(item) for item in obj]
        else:
            return obj

    results_serializable = convert_types(results)

    output_file = results_dir / "quaternion_bridge_2i_audit.json"
    with open(output_file, "w") as f:
        json.dump(results_serializable, f, indent=2)

    print(f"Results saved to: {output_file}")
    print()

    # Final summary
    print("=" * 70)
    print("AUDIT SUMMARY")
    print("=" * 70)
    print(f"Vierbein construction: {'✓ PASS' if vierbein_metadata['positive_definite_count'] == 64 else '⚠ PARTIAL'}")
    print(f"SU(2) connection: {'✓ PASS' if len(su2_connection) == 192 else '✗ FAIL'}")
    print(f"2I group: {'✓ PASS' if elements_2I.shape[0] == 120 else '✗ FAIL'}")
    print(f"Link 2I proximity: {mean_dist:.6f}")
    print(f"Holonomy 2I proximity: {holonomy_mean:.6f}")
    print(f"Max curvature: {np.degrees(ca['max']):.2f}°")
    print()

    if mean_dist < 0.1:
        print("✓ CRITICAL FINDING: Strong spontaneous alignment with 2I")
        print("  This supports the hypothesis that the geometric connection")
        print("  naturally approaches the binary icosaededral structure.")
    elif mean_dist < 0.5:
        print("⚠ MODERATE FINDING: Partial alignment with 2I")
        print("  Further analysis needed to determine significance.")
    else:
        print("✗ NEGATIVE FINDING: No spontaneous 2I alignment")
        print("  The geometric connection does not approach 2I.")
        print("  This challenges the hypothesis of automatic 2I emergence.")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
