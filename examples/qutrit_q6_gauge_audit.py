import json
import numpy as np
from src.qutrit_q6_coupling import build_q6_qutrit_chain, q6_qutrit_gauge_covariance_residual, q6_qutrit_spectrum

def diagonal_su3(a, b):
    return np.diag([np.exp(1j*a), np.exp(1j*b), np.exp(-1j*(a+b))])

if __name__ == '__main__':
    chain = build_q6_qutrit_chain(4, j_a=1.0, j_b=0.37)
    gauge = {v: diagonal_su3(0.07*(i+1), -0.11*(i+1)) for i, v in enumerate(chain.site_vertices)}
    print(json.dumps({'site_vertices': list(chain.site_vertices), 'fibonacci_word': chain.word, 'gauge_covariance_residual': q6_qutrit_gauge_covariance_residual(chain, gauge), 'spectrum_and_entanglement': q6_qutrit_spectrum(chain)}, indent=2))
