import numpy as np


def compute_fine_structure_derivation():
  # Invarianti geometrici puri già estratti dal nostro operatore
  W_geo = 0.188562  # Quanto d'azione medio per faccia ibrida
  hbar_geo = 2.828427  # 2 * sqrt(2)
  n_edges = 12
  n_faces = 24
  vol_fraction = 2.0 / 3.0  # Volume occupato dalla stella octangula nel cubo

  # Valore di riferimento reale di alpha^-1
  alpha_inv_real = 137.035999

  # Tentativo 1: Accoppiamento basato sul rapporto fra spigoli, facce e azione
  # Usando la chiusura geometrica del gruppo di simmetria (es. volume del cubo = 8)
  # alpha^-1 ~ (4 * pi^3 / vol_fraction) corretto per l'azione d'interfaccia?
  # Proviamo la relazione topologica pura con la chiralità:
  
  # Relazione geometrica nativa basata sui 12 spigoli e la radice di 2:
  # L'interazione tra i 2 tetraedri in un cubo di volume 8
  term_base = (8.0 / vol_fraction) * (n_faces / n_edges) # 8 / (2/3) * 2 = 24 * 3 = 72
  
  # Incorporando l'azione geometrica W_geo e la metrica di Planck geometrica hbar_geo:
  # Cerchiamo l'espressione adimensionale esatta:
  alpha_inv_calc = (16.0 * (np.pi**2) * hbar_geo) / (W_geo * n_edges)

  print("=== DERIVAZIONE GEOMETRICA DI ALPHA ===")
  print(f"Valore reale di α^(-1):          {alpha_inv_real:.6f}")
  print(f"Valore derivato α^(-1) (geo):    {alpha_inv_calc:.6f}")
  
  # Differenza percentuale
  diff_pct = abs(alpha_inv_calc - alpha_inv_real) / alpha_inv_real * 100
  print(f"Scostamento percentuale:         {diff_pct:.4f}%")

  # Verifichiamo un'altra combinazione basata sui pesi di Schmidt (1/4 e 1/12)
  # 1/4 + 1/12 = 1/3
  schmidt_sum = (1.0 / 4.0) + (1.0 / 12.0) # 1/3
  alpha_inv_schmidt = (4.0 * np.pi**3) / (schmidt_sum * W_geo)
  print(f"Valore alternativo (Schmidt):    {alpha_inv_schmidt:.6f}")

if __name__ == "__main__":
  compute_fine_structure_derivation()
