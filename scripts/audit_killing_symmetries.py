"""Complete audit of RGS vierbein construction and Killing symmetries.

This script:
1. Loads RGS matrices from NPZ archive
2. Constructs vierbein from covariance matrices
3. Reconstructs metric
4. Searches for Killing vector fields
5. Verifies metric symmetries
6. Analyzes the symmetry algebra
7. Generates comprehensive audit report
"""

import sys
import numpy as np
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent))

from vierbein.audit import (
    audit_vierbein,
    audit_metric,
    audit_killing_symmetries,
    print_audit_report,
)
from load_rgs_matrices import load_rgs_matrices


def main():
    """Main execution function."""
    print("=" * 70)
    print("RGS-PSU KILLING SYMMETRY AUDIT")
    print("=" * 70)
    print()

    # Load RGS matrices
    print("Loading RGS matrices...")
    matrices, metadata = load_rgs_matrices()
    print(f"Loaded {len(matrices)} arrays from {metadata['file_path']}")
    print(f"SHA-256: {metadata['sha256']}")
    print()

    # Extract required data
    if "X_Oh_reconstructed" not in matrices:
        print("ERROR: X_Oh_reconstructed not found in matrices")
        sys.exit(1)

    X = matrices["X_Oh_reconstructed"]
    print(f"Vertex positions X: shape {X.shape}")
    print()

    # Historical parameters (to be verified)
    ell_squared = 7.035514834154695
    c0 = 43.96495130838031

    print("Using historical parameters:")
    print(f"  ℓ² = {ell_squared}")
    print(f"  c₀ = {c0}")
    print()

    # Step 1: Vierbein construction audit
    print("STEP 1: VIERBEIN CONSTRUCTION")
    print("-" * 70)
    vierbein_audit = audit_vierbein(X, ell_squared, c0, tol=1e-12)

    vc = vierbein_audit["vierbein_construction"]
    print(f"Positive definite matrices: {vc['positive_definite_count']}/{vc['total_matrices']}")

    if vc["positive_definite_count"] < vc["total_matrices"]:
        print("WARNING: Not all covariance matrices are positive definite!")
        print("This may affect vierbein construction.")
    print()

    # Check metric identity
    mi = vierbein_audit["metric_identity"]
    print(f"Metric identity g = C/c₀:")
    print(f"  Max residual: {mi['max_residual']:.2e}")
    print(f"  Identity holds: {mi['identity_holds']}")

    if not mi["identity_holds"]:
        print("WARNING: Metric identity does not hold within tolerance!")
    print()

    # Step 2: Reconstruct metric for further analysis
    print("STEP 2: METRIC RECONSTRUCTION")
    print("-" * 70)
    from vierbein.construction import build_vierbein, reconstruct_metric

    E, _ = build_vierbein(X, ell_squared, c0)
    g = reconstruct_metric(E)

    print(f"Vierbein shape: {E.shape}")
    print(f"Metric shape: {g.shape}")
    print()

    # Step 3: Metric properties audit
    print("STEP 3: METRIC PROPERTIES AUDIT")
    print("-" * 70)
    metric_audit = audit_metric(g, X, tol=1e-10)

    mp = metric_audit["symmetry"]
    print(f"Symmetry check:")
    print(f"  Max error: {mp['max_error']:.2e}")
    print(f"  All symmetric: {mp['all_symmetric']}")

    pd = metric_audit["positive_definiteness"]
    print(f"Positive definiteness:")
    print(f"  Count: {pd['count']}/{pd['total']}")
    print(f"  Min eigenvalue: {pd['min_eigenvalue']:.2e}")
    print(f"  All positive: {pd['all_positive']}")
    print()

    # Step 4: Killing symmetry audit
    print("STEP 4: KILLING SYMMETRY ANALYSIS")
    print("-" * 70)

    # Use Q6 adjacency if available
    adjacency = matrices.get("Q6_adjacency", None)
    if adjacency is not None:
        print(f"Using Q6 adjacency matrix: shape {adjacency.shape}")
    else:
        print("No adjacency matrix provided, will compute from distances")
    print()

    killing_audit = audit_killing_symmetries(X, g, adjacency, tol=1e-8)

    print(f"Killing vectors found: {killing_audit['num_killing_vectors']}")
    print(f"Expected (Euclidean 3D): {killing_audit['expected_dimension']}")
    print(f"Algebra closed: {killing_audit['algebra_closed']}")

    if killing_audit["killing_residuals"]:
        print(f"Killing equation residuals:")
        for i, res in enumerate(killing_audit["killing_residuals"]):
            print(f"  Field {i}: {res:.2e}")
    print()

    # Step 5: Analysis of results
    print("STEP 5: INTERPRETATION")
    print("-" * 70)

    num_killing = killing_audit["num_killing_vectors"]

    if num_killing == 6:
        print("✓ Full Euclidean symmetry group (3 translations + 3 rotations)")
        print("  The reconstructed metric has the expected continuous symmetries.")
    elif num_killing >= 3:
        print(f"⚠ Partial symmetry: {num_killing} Killing vectors found")
        print("  This suggests some but not all Euclidean symmetries are present.")
    elif num_killing > 0:
        print(f"⚠ Limited symmetry: {num_killing} Killing vector(s) found")
        print("  The metric has restricted continuous symmetries.")
    else:
        print("✗ No exact Killing vectors found within tolerance")
        print("  The reconstructed metric lacks continuous symmetries.")
        print("  This may indicate:")
        print("    - Discrete-only symmetries (Q6, 2I)")
        print("    - Approximate symmetries (increase tolerance)")
        print("    - Genuine lack of continuous isometries")
    print()

    # Step 6: Connection to discrete symmetries
    print("STEP 6: CONNECTION TO DISCRETE SYMMETRIES")
    print("-" * 70)
    print("The discrete symmetries available are:")
    print("  - Q6 automorphisms (combinatorial symmetries of hypercube)")
    print("  - Binary icosahedral group 2I (120 quaternionic rotations)")
    print("  - O_h symmetry (cubic symmetry of X_Oh_reconstructed)")
    print()
    print("Key question: Do these discrete symmetries extend to continuous")
    print("Killing vectors of the reconstructed metric?")
    print()

    if num_killing == 6:
        print("ANSWER: YES - Full Euclidean symmetry present")
        print("  → Discrete symmetries are subgroups of continuous symmetry")
    elif num_killing > 0:
        print("ANSWER: PARTIAL - Some continuous symmetries present")
        print("  → Further analysis needed to map discrete → continuous")
    else:
        print("ANSWER: NO - No continuous symmetries found")
        print("  → Discrete symmetries do NOT extend to Killing vectors")
        print("  → This is a significant scientific result to report")
    print()

    # Step 7: Save results
    print("STEP 7: SAVING RESULTS")
    print("-" * 70)

    results_dir = Path("results/json")
    results_dir.mkdir(parents=True, exist_ok=True)

    import json

    results = {
        "metadata": {
            "file_hash": metadata["sha256"],
            "ell_squared": ell_squared,
            "c0": c0,
            "num_vertices": int(X.shape[0]),
        },
        "vierbein_audit": vierbein_audit,
        "metric_audit": metric_audit,
        "killing_audit": killing_audit,
        "interpretation": {
            "num_killing_vectors": num_killing,
            "has_full_euclidean_symmetry": num_killing == 6,
            "has_partial_symmetry": 0 < num_killing < 6,
            "has_no_symmetry": num_killing == 0,
        },
    }

    # Convert numpy types to Python types for JSON serialization
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

    output_file = results_dir / "killing_symmetry_audit.json"
    with open(output_file, "w") as f:
        json.dump(results_serializable, f, indent=2)

    print(f"Results saved to: {output_file}")
    print()

    # Final summary
    print("=" * 70)
    print("AUDIT SUMMARY")
    print("=" * 70)
    print(f"Vierbein construction: {'✓ PASS' if vc['positive_definite_count'] == vc['total_matrices'] else '⚠ PARTIAL'}")
    print(f"Metric identity: {'✓ PASS' if mi['identity_holds'] else '✗ FAIL'}")
    print(f"Metric positive definite: {'✓ PASS' if pd['all_positive'] else '✗ FAIL'}")
    print(f"Killing vectors: {num_killing}/6 expected")
    print()

    if num_killing == 6:
        print("✓ AUDIT PASSED: Full continuous symmetry confirmed")
        sys.exit(0)
    elif num_killing > 0:
        print("⚠ AUDIT WARNING: Partial symmetry detected")
        sys.exit(0)
    else:
        print("✗ AUDIT FAILED: No continuous symmetries found")
        print("  This is a scientifically significant result.")
        sys.exit(0)


if __name__ == "__main__":
    main()
