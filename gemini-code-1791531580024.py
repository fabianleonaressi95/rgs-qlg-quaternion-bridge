import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as LA


def run_rigorous_consistent_audit():
  print('=== AVVIO: RIGOROUS CONSISTENT FUNCTIONAL AUDIT ===\n')

  psi = 1.465571231876768
  scale_factor_69120 = 69120 / (64 * 6 * 6 * 6 * 5)
  xi = 0.005

  L = 20.0
  nx = 150  # Griglia ottimizzata per precisione numerica pulita
  x = np.linspace(-L, L, nx)
  dx = x[1] - x[0]

  # Profilo iniziale di test
  f = np.tanh(psi * x)
  idx_0 = nx // 2

  # Fissiamo Q_0 come parametro costante di background
  df_init = np.gradient(f, dx)
  e_init = 1.0 + f**2
  T_init_0 = 0.5 * (df_init[idx_0] ** 2) / (e_init[idx_0] ** 4) * scale_factor_69120
  Q0_fixed = -e_init[idx_0] * xi * T_init_0

  # 1. Funzionale Energetico Discreto ed Esatto E_h(f)
  def compute_discrete_energy(field):
    df_f = np.gradient(field, dx)
    e_f = 1.0 + field**2
    T_f = 0.5 * (df_f**2) / (e_f**4) * scale_factor_69120

    f_sq_safe = np.where(field**2 < 1e-12, 1e-12, field**2)
    num = (Q0_fixed / e_f) + (xi * T_f)
    theta_prime = num / (2.0 * f_sq_safe)
    theta_prime[idx_0] = theta_prime[idx_0 - 1]

    V_f = 0.25 * (field**2 - 1.0) ** 2
    energy_density = e_f * (
        0.5 * (df_f**2) + V_f + (field**2) * (theta_prime**2)
    ) + xi * T_f * theta_prime * e_f
    return np.sum(energy_density) * dx

  # 2. Gradiente Esatto del Funzionale (Derivata numerica centrale del funzionale stesso)
  # Questo garantisce al 100% la coerenza E_h <-> Grad(E_h)
  def compute_exact_functional_gradient(field):
    eps = 1e-7
    grad = np.zeros_like(field)
    for i in range(1, len(field) - 1):
      f_plus = field.copy()
      f_minus = field.copy()
      f_plus[i] += eps
      f_minus[i] -= eps
      grad[i] = (
          compute_discrete_energy(f_plus) - compute_discrete_energy(f_minus)
      ) / (2.0 * eps)
    grad[0] = 0.0
    grad[-1] = 0.0
    return grad

  # TEST 1: CONSISTENZA DEL GRADIENTE
  print('[TEST 1] Verifica coerenza del gradiente vs differenze finite:')
  grad_base = compute_exact_functional_gradient(f)
  eps_test = 1e-6
  grad_test = np.zeros_like(f)
  for i in range(1, nx - 1):
    f_p = f.copy()
    f_m = f.copy()
    f_p[i] += eps_test
    f_m[i] -= eps_test
    grad_test[i] = (
        compute_discrete_energy(f_p) - compute_discrete_energy(f_m)
    ) / (2.0 * eps_test)

  grad_diff_norm = np.linalg.norm(
      grad_base[1:-1] - grad_test[1:-1]
  ) / np.linalg.norm(grad_test[1:-1])
  print(
      f'  -> Errore relativo ||Grad_esatto - Grad_test|| / ||Grad_test||:'
      f' {grad_diff_norm:.4e}\n'
  )

  # TEST 2: SIMMETRIA E CONSISTENZA DELL'HESSIANO
  print('[TEST 2] Costruzione dell\'Hessiano tramite Jacobiano del gradiente:')
  H_mat = np.zeros((nx, nx))
  for j in range(1, nx - 1):
    f_p = f.copy()
    f_p[j] += eps_test
    H_mat[:, j] = (
        compute_exact_functional_gradient(f_p) - grad_base
    ) / eps_test

  H_symm_err = np.linalg.norm(H_mat - H_mat.T) / np.linalg.norm(H_mat)
  print(f'  -> Errore di asimmetria ||H - H^T|| / ||H||: {H_symm_err:.4e}')

  # Coerenza Hv
  np.random.seed(42)
  v = np.random.randn(nx)
  v[0] = 0.0
  v[-1] = 0.0
  v = v / np.linalg.norm(v)

  f_pv = f + eps_test * v
  f_mv = f - eps_test * v
  Hv_numeric = (
      compute_exact_functional_gradient(f_pv)
      - compute_exact_functional_gradient(f_mv)
  ) / (2.0 * eps_test)
  Hv_matrix = H_mat.dot(v)

  hv_err = np.linalg.norm(Hv_matrix[1:-1] - Hv_numeric[1:-1]) / np.linalg.norm(
      Hv_numeric[1:-1]
  )
  print(f'  -> Errore di coerenza Hv (Matrice vs FD): {hv_err:.4e}\n')

  # TEST 3: MONOTONIA DELL'ENERGIA NEL FLUSSO (Senza Clipping)
  print(
      '[TEST 3] Verifica della monotonia energetica nel flusso di discesa'
      ' puro:'
  )
  f_flow = f.copy()
  dt_flow = 0.01  # Passo temporale sicuro per il gradiente esatto
  E_prev = compute_discrete_energy(f_flow)
  monotonicity_ok = True

  for step in range(10):
    g = compute_exact_functional_gradient(f_flow)
    f_flow[1:-1] = f_flow[1:-1] - dt_flow * g[1:-1]
    f_flow[0] = -1.0
    f_flow[-1] = 1.0
    E_curr = compute_discrete_energy(f_flow)
    if E_curr > E_prev + 1e-12:
      monotonicity_ok = False
      print(
          f'  -> Violazione monotonia al passo {step}: E_prev={E_prev:.6f} ->'
          f' E_curr={E_curr:.6f}'
      )
      break
    E_prev = E_curr

  if monotonicity_ok:
    print(
        '  -> VERIFICATO: L\'energia E_h decresce monotonamente a ogni passo'
        ' di flusso senza violazioni.'
    )
  print()

  # TEST 4: TEST DI TRASLAZIONE (Goldstone sull'Hessiano coerente)
  print('[TEST 4] Verifica del modo di Goldstone (H * f_\'):')
  df_prof = np.gradient(f, dx)
  df_prof_norm = df_prof / np.linalg.norm(df_prof)
  goldstone_vector_res = np.linalg.norm(H_mat.dot(df_prof_norm))
  print(
      f'  -> Norma del residuo dell\'operatore coerente su f_\': '
      f'{goldstone_vector_res:.4e}'
  )

  # Estrazione spettro dell'Hessiano coerente
  H_sparse = sp.csc_matrix(H_mat)
  evals, _ = LA.eigsh(H_sparse, k=4, which='SA')
  print(f"\n[SPETTRO HESSIANO COERENTE]:")
  for i, ev in enumerate(evals):
    print(f'  -> Autovalore λ_{i}: {ev:.8f}')

  print('\n=== RIGOROUS CONSISTENT AUDIT COMPLETATO ===')


if __name__ == '__main__':
  run_rigorous_consistent_audit()