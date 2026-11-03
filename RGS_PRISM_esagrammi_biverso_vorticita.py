# ================================================================
# RGS–PRISM | 64 ESAGRAMMI -> BIVERSO -> AUDIT VORTICITÀ DISCRETA
#
# Obiettivo:
#   A) costruire il grafo Q6 = Z2^6;
#   B) costruire il biverso Q6 x Z2 e lo scambio dei fogli;
#   C) calcolare olonomie su quadrati elementari se è disponibile
#      una connessione U_i(g,s) derivata dai dati geometrici;
#   D) se la connessione non esiste, NON inventare curvature:
#      mostrare che la topologia da sola non determina la vorticità.
#
# Convenzione:
#   esagramma g = intero 0..63, bit i = linea i;
#   transizione i: g -> g XOR (1 << i);
#   foglio s = 0,1.
#
# File opzionali:
#   CONNECTION_NPY: array complesso/float di shape (2,64,6),
#     U[s,g,i] = fase di link orientato da g a g XOR 2**i.
#     Può contenere valori complessi unitari, oppure fasi reali A
#     in radianti (se reali, lo script converte exp(1j*A)).
#   In assenza del file, la connessione resta NON SPECIFICATA.
#
# Nota: una connessione costante U=1 darebbe curvatura nulla.
# Non è una scoperta fisica; è soltanto il caso piatto di controllo.
# ================================================================

import os
import json
import numpy as np

OUT = "/kaggle/working/RGS_PRISM_ESAGRAMMI_VORTICITA"
os.makedirs(OUT, exist_ok=True)

# Modifica questo percorso SOLO se possiedi una connessione
# geometricamente derivata e salvata in un file.
CONNECTION_NPY = None
# Esempio:
# CONNECTION_NPY = "/kaggle/input/.../U_esagrammi_biverso.npy"

N_G = 64
N_LINE = 6
N_SHEETS = 2
N_BI = N_G * N_SHEETS
TOL = 1e-10

def bit_flip(g, i):
    return int(g) ^ (1 << int(i))

def node_id(g, s):
    return int(s) * N_G + int(g)

# ------------------------------------------------
# 1. Grafo Q6 e audit combinatorio
# ------------------------------------------------
edges = []
adj = [[] for _ in range(N_G)]
for g in range(N_G):
    for i in range(N_LINE):
        h = bit_flip(g, i)
        adj[g].append(h)
        if g < h:
            edges.append((g, h, i))

degree = np.array([len(x) for x in adj])
assert len(edges) == 192
assert np.all(degree == 6)

# Ogni coppia di coordinate i<j genera 16*? quadrati distinti:
# per ogni asse i,j, i quadrati sono indicizzati dai bit restanti:
# 2^(6-2)=16; totale C(6,2)*16 = 240 quadrati.
squares = []
for i in range(N_LINE):
    for j in range(i + 1, N_LINE):
        for g in range(N_G):
            # Canonical representative: i e j bit a zero.
            if (g & (1 << i)) or (g & (1 << j)):
                continue
            gi = bit_flip(g, i)
            gj = bit_flip(g, j)
            gij = bit_flip(gi, j)
            squares.append((g, gi, gij, gj, i, j))

assert len(squares) == 240

# ------------------------------------------------
# 2. Biverso e scambio dei fogli
# ------------------------------------------------
# S(g,0)=(g,1), S(g,1)=(g,0)
swap = np.array([node_id(g, 1-s) for s in range(2) for g in range(64)])
assert np.array_equal(swap[swap], np.arange(N_BI))

# Proiettori sullo spazio vettoriale delle funzioni sui 128 nodi.
# Non materializziamo matrici 128x128: testiamo l'azione direttamente.
def apply_swap(v):
    v = np.asarray(v).reshape(2, N_G, -1)
    return v[::-1].copy().reshape(2*N_G, -1)

rng = np.random.default_rng(20261021)
v = rng.normal(size=(N_BI, 3))
Sv = apply_swap(v)
S2v = apply_swap(Sv)
swap_involution_res = float(np.linalg.norm(S2v-v))

vplus = (v + Sv) / 2
vminus = (v - Sv) / 2
decomp_res = float(np.linalg.norm(vplus + vminus - v))
orthogonality = float(np.linalg.norm(vplus.T @ vminus))
plus_swap_res = float(np.linalg.norm(apply_swap(vplus)-vplus))
minus_swap_res = float(np.linalg.norm(apply_swap(vminus)+vminus))

# ------------------------------------------------
# 3. Caricamento connessione, se esiste
# ------------------------------------------------
conn_status = "NOT_SPECIFIED"
U = None
curvature_results = None

if CONNECTION_NPY is not None and os.path.exists(CONNECTION_NPY):
    raw = np.load(CONNECTION_NPY)
    if raw.shape != (2, 64, 6):
        raise ValueError(
            f"Connessione shape {raw.shape}; attesa (2,64,6)."
        )

    if np.iscomplexobj(raw):
        U = raw.astype(np.complex128)
    else:
        # Una matrice reale è interpretata come fase A in radianti.
        U = np.exp(1j * raw.astype(np.float64))

    # Per ogni link, U deve essere una fase unitaria.
    unit_modulus_res = float(np.max(np.abs(np.abs(U)-1)))
    if unit_modulus_res > 1e-7:
        raise ValueError(
            "Connessione complessa non unitaria: max ||U|-1| = "
            f"{unit_modulus_res}"
        )

    # Il link inverso è il coniugato del link diretto.
    # Calcoliamo il Wilson loop quadrato:
    # U_i(g) U_j(g^i) U_i(g^j)^* U_j(g)^*
    # in ogni foglio.
    W = np.empty((2, len(squares)), dtype=np.complex128)
    for s in range(2):
        for k, (g, gi, gij, gj, i, j) in enumerate(squares):
            W[s, k] = (
                U[s, g, i]
                * U[s, gi, j]
                * np.conj(U[s, gj, i])
                * np.conj(U[s, g, j])
            )

    phases = np.angle(W)
    curvature_abs = np.abs(phases)
    flat_counts = np.sum(curvature_abs < TOL, axis=1)

    # Test di curvatura opposta tra fogli: W_1 ~= conj(W_0).
    opposite_res = float(
        np.max(np.abs(W[1] - np.conj(W[0])))
    )
    same_res = float(np.max(np.abs(W[1] - W[0])))

    curvature_results = {
        "connection_path": CONNECTION_NPY,
        "unit_modulus_max_residual": unit_modulus_res,
        "plaquettes_per_sheet": int(len(squares)),
        "flat_plaquettes_sheet0": int(flat_counts[0]),
        "flat_plaquettes_sheet1": int(flat_counts[1]),
        "nonflat_plaquettes_sheet0": int(len(squares)-flat_counts[0]),
        "nonflat_plaquettes_sheet1": int(len(squares)-flat_counts[1]),
        "max_abs_plaquette_phase_sheet0": float(np.max(curvature_abs[0])),
        "max_abs_plaquette_phase_sheet1": float(np.max(curvature_abs[1])),
        "rms_plaquette_phase_sheet0": float(np.sqrt(np.mean(phases[0]**2))),
        "rms_plaquette_phase_sheet1": float(np.sqrt(np.mean(phases[1]**2))),
        "opposite_sheet_curvature_residual": opposite_res,
        "same_sheet_curvature_residual": same_res,
        "opposite_sheet_curvature_within_1e-8": bool(opposite_res < 1e-8),
        "same_sheet_curvature_within_1e-8": bool(same_res < 1e-8),
        "interpretation": (
            "Curvatura calcolata dalla connessione fornita. "
            "La sua origine geometrica deve essere documentata separatamente."
        )
    }
    conn_status = "LOADED_AND_AUDITED"
else:
    # Non inseriamo U=exp(iA) artificiale per dichiarare vortici.
    conn_status = "NOT_SPECIFIED_NO_VORTICITY_CLAIM"

# ------------------------------------------------
# 4. Audit dimensionale e interpretazione
# ------------------------------------------------
report = {
    "model": "Q6 esagrams x Z2 two-sheet biverso",
    "dimensions": {
        "esagrams": N_G,
        "lines_per_esagram": N_LINE,
        "base_edges_undirected": len(edges),
        "elementary_square_plaquettes_per_sheet": len(squares),
        "sheets": N_SHEETS,
        "biverso_states": N_BI,
        "symmetric_sector_dimension": N_G,
        "antisymmetric_sector_dimension": N_G
    },
    "combinatorial_checks": {
        "degree_min": int(degree.min()),
        "degree_max": int(degree.max()),
        "swap_square_identity_residual": swap_involution_res,
        "Pplus_plus_Pminus_identity_residual": decomp_res,
        "Pplus_Pminus_inner_product_norm": orthogonality,
        "S_Pplus_equals_Pplus_residual": plus_swap_res,
        "S_Pminus_equals_minus_Pminus_residual": minus_swap_res
    },
    "connection": {
        "status": conn_status,
        "curvature": curvature_results
    },
    "conclusions": {
        "Q6_biverso_combinatorics": "VERIFIED_NUMERICALLY",
        "vorticity_from_esagrams_alone": "NOT_DETERMINED",
        "double_toroidal_physical_interpretation": "OPEN_HYPOTHESIS",
        "required_next_input": (
            "Una connessione orientata U[s,g,i] derivata indipendentemente "
            "dalla geometria RGS (o una regola geometrica equivalente)."
        )
    }
}

np.save(os.path.join(OUT, "Q6_edges.npy"), np.array(edges, dtype=np.int64))
np.save(os.path.join(OUT, "Q6_squares.npy"), np.array(squares, dtype=np.int64))
np.save(os.path.join(OUT, "biverso_swap_index.npy"), swap)

with open(os.path.join(OUT, "Q6_biverso_vorticity_audit.json"), "w") as f:
    json.dump(report, f, indent=2)

print("=" * 72)
print("RGS–PRISM | ESAGRAMMI / BIVERSO / VORTICITÀ")
print("=" * 72)
print(json.dumps(report, indent=2))
print("\nOutput:", OUT)
