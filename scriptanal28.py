import numpy as np


def refined_alpha_derivation():
  W_geo = 0.188562  # Quanto d'azione medio per faccia ibrida
  n_edges = 12
  n_faces = 24
  vol_fraction = 2.0 / 3.0  # 2/3
  alpha_inv_real = 137.035999

  # Correzione basata sul bilancio di simmetria tra il volume del cubo (8),
  # i 12 spigoli e l'azione di interfaccia W_geo:
  # Cerchiamo il fattore scalare puro basato sulla geometria solida della stella.
  
  # Tentativo di calcolo pulito con la metrica d'interfaccia pura:
  # alpha^-1 = ( Volume Cubo * n_edges^2 ) / ( 4 * pi * W_geo * vol_fraction )? No, testiamo la proporzione esatta:
  
  # Relazione geometrica bilanciata sui 12 spigoli e le 24 facce ibride:
  # 137 è vicinissimo a (12^2) + qualcosa, oppure 4 * pi * 11 ... 
  # Verifichiamo la formula di risonanza geometrica pura:
  
  val_test = (4.0 * np.pi * (n_edges ** 2)) / (W_geo * vol_fraction * np.sqrt(2))
  
  print("=== RIFINITURA GEOMETRICA DI ALPHA ===")
  print(f"Valore reale di α^(-1):        {alpha_inv_real:.6f}")
  print(f"Valore test geometrico:        {val_test:.6f}")
  
  # Troviamo esattamente quale fattore numerico K chiude l'equazione a 137.036:
  # K = alpha_inv_real * (W_geo * vol_fraction * np.sqrt(2)) / (4 * pi * 144)
  K_factor = (alpha_inv_real * W_geo * vol_fraction * np.sqrt(2)) / (4.0 * np.pi * (n_edges ** 2))
  print(f"Fattore di chiusura K emerso:  {K_factor:.6f}")

if __name__ == "__main__":
  refined_alpha_derivation()
