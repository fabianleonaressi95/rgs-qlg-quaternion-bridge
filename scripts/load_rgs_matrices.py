"""Load and validate RGS-PSU matrices from NPZ archive.

This script loads the RGS_PSU_all_matrices.npz file and provides
access to the historical matrices with validation checks.
"""

import numpy as np
from pathlib import Path
from typing import Dict, Tuple
import hashlib


def compute_sha256(filepath: Path) -> str:
    """Compute SHA-256 hash of a file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def load_rgs_matrices(
    data_path: Path = Path("data/rgs_188d/RGS_PSU_all_matrices.npz"),
) -> Tuple[Dict[str, np.ndarray], Dict[str, str]]:
    """Load RGS matrices from NPZ archive with validation.

    Returns:
        matrices: Dictionary mapping keys to numpy arrays
        metadata: Dictionary with file metadata (hash, shape, dtype)
    """
    if not data_path.exists():
        raise FileNotFoundError(f"RGS matrices not found at {data_path}")

    # Compute file hash for provenance
    file_hash = compute_sha256(data_path)

    # Load all arrays
    data = np.load(data_path)

    matrices = {}
    metadata = {
        "file_path": str(data_path),
        "sha256": file_hash,
        "keys": list(data.keys()),
    }

    for key in data.keys():
        arr = data[key]
        matrices[key] = arr
        metadata[f"{key}_shape"] = str(arr.shape)
        metadata[f"{key}_dtype"] = str(arr.dtype)
        metadata[f"{key}_size_mb"] = arr.nbytes / (1024 * 1024)

    return matrices, metadata


def validate_rgs_matrices(matrices: Dict[str, np.ndarray]) -> Dict[str, bool]:
    """Run validation checks on RGS matrices.

    Returns:
        validation_results: Dictionary mapping check names to pass/fail
    """
    results = {}

    # Check X_Oh_reconstructed
    if "X_Oh_reconstructed" in matrices:
        X = matrices["X_Oh_reconstructed"]
        results["X_shape_is_64x3"] = X.shape == (64, 3)
        results["X_is_finite"] = np.all(np.isfinite(X))
        results["X_dtype_is_float64"] = X.dtype == np.float64

    # Check Q_tangent
    if "Q_tangent" in matrices:
        Q = matrices["Q_tangent"]
        results["Q_shape_is_192x188"] = Q.shape == (192, 188)
        results["Q_is_finite"] = np.all(np.isfinite(Q))
        results["Q_dtype_is_float64"] = Q.dtype == np.float64

    # Check Q6 structures
    if "Q6_states" in matrices:
        states = matrices["Q6_states"]
        results["Q6_states_shape_is_64x6"] = states.shape == (64, 6)
        results["Q6_states_binary"] = np.all(
            (states == 0) | (states == 1)
        )

    if "Q6_adjacency" in matrices:
        adj = matrices["Q6_adjacency"]
        results["Q6_adjacency_shape_is_64x64"] = adj.shape == (64, 64)
        results["Q6_adjacency_symmetric"] = np.allclose(adj, adj.T)
        results["Q6_adjacency_degree_6"] = np.allclose(
            np.sum(adj > 0, axis=1), 6
        )

    if "Q6_laplacian" in matrices:
        lap = matrices["Q6_laplacian"]
        results["Q6_laplacian_shape_is_64x64"] = lap.shape == (64, 64)
        results["Q6_laplacian_symmetric"] = np.allclose(lap, lap.T)

    # Check operator matrices
    if "K0" in matrices:
        K0 = matrices["K0"]
        results["K0_shape_is_188x188"] = K0.shape == (188, 188)
        results["K0_symmetric"] = np.allclose(K0, K0.T)
        results["K0_finite"] = np.all(np.isfinite(K0))

    if "P_G" in matrices:
        P_G = matrices["P_G"]
        results["P_G_shape_is_188x188"] = P_G.shape == (188, 188)
        # Check if it's a projector: P^2 = P
        P_G_squared = P_G @ P_G
        results["P_G_is_projector"] = np.allclose(P_G, P_G_squared, atol=1e-10)
        results["P_G_symmetric"] = np.allclose(P_G, P_G.T)

    return results


def print_validation_report(
    matrices: Dict[str, np.ndarray],
    metadata: Dict[str, str],
    validation: Dict[str, bool],
):
    """Print a formatted validation report."""
    print("=" * 70)
    print("RGS-PSU MATRICES VALIDATION REPORT")
    print("=" * 70)
    print()

    print("FILE METADATA")
    print("-" * 70)
    print(f"Path: {metadata['file_path']}")
    print(f"SHA-256: {metadata['sha256']}")
    print(f"Keys: {', '.join(metadata['keys'])}")
    print()

    print("ARRAY SUMMARY")
    print("-" * 70)
    for key in matrices.keys():
        shape = metadata[f"{key}_shape"]
        dtype = metadata[f"{key}_dtype"]
        size_mb = metadata[f"{key}_size_mb"]
        print(f"{key:30s} {shape:15s} {dtype:10s} {size_mb:8.2f} MB")
    print()

    print("VALIDATION CHECKS")
    print("-" * 70)
    passed = 0
    failed = 0
    for check, result in validation.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status:8s} {check}")
        if result:
            passed += 1
        else:
            failed += 1

    print()
    print(f"Total: {passed + failed} checks")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print("=" * 70)


def main():
    """Main entry point for script execution."""
    import sys

    # Load matrices
    print("Loading RGS matrices...")
    matrices, metadata = load_rgs_matrices()
    print(f"Loaded {len(matrices)} arrays")
    print()

    # Validate
    print("Running validation checks...")
    validation = validate_rgs_matrices(matrices)
    print()

    # Print report
    print_validation_report(matrices, metadata, validation)

    # Return exit code
    if all(validation.values()):
        print("\n✓ All validation checks passed")
        sys.exit(0)
    else:
        print("\n✗ Some validation checks failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
