def calculate_rigorous_diagnostics(H, coordinates, filling_fraction=0.7):
    """
    Calcola il marcatore locale di Chern e verifica la stabilità del gap.
    """
    N = H.shape[0]
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    evals, eigenvectors = LA.eigh(H_dense)
    
    num_occupied = int(N * filling_fraction)
    if num_occupied <= 0 or num_occupied >= N:
        return 0.0, 0.0, False
        
    mean_level_spacing = np.mean(np.diff(evals))
    # CORRETTO: eval -> evals
    local_gap = evals[num_occupied] - evals[num_occupied - 1] if 0 < num_occupied < N else 0.0
    is_true_gap = local_gap > (1.2 * mean_level_spacing)
    
    V_occ = eigenvectors[:, :num_occupied]
    P = V_occ @ V_occ.conj().T
    
    min_coords = coordinates.min(axis=0)
    max_coords = coordinates.max(axis=0)
    scaled_coords = (coordinates - min_coords) / (max_coords - min_coords + 1e-12)
    
    X = np.diag(scaled_coords[:, 0])
    Y = np.diag(scaled_coords[:, 1])
    
    P_X_P_Y_P = P @ X @ P @ Y @ P
    local_chern = -4.0 * np.pi * np.imag(np.diagonal(P_X_P_Y_P))
    
    center = coordinates.mean(axis=0)
    distances = np.linalg.norm(coordinates - center, axis=1)
    bulk_mask = distances < (0.35 * distances.max())
    
    mean_bulk_chern = np.mean(local_chern[bulk_mask]) if np.any(bulk_mask) else 0.0
    
    return mean_bulk_chern, local_gap, is_true_gap
