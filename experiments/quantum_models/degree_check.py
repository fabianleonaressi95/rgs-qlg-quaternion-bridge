import numpy as np

def T_pell(a, b):
    return ((2 * a + b) % 8, a % 8)

def F(n):
    a = n % 8
    b = n // 8
    a2, b2 = T_pell(a, b)
    return a2 + 8 * b2

perm = [F(n) for n in range(64)]
K = np.zeros((64, 64))
for x in range(64):
    K[x, perm[x]] = 1.0

def bits(n):
    return [(n >> i) & 1 for i in range(6)]

W_raw = np.zeros((64, 64))
for x in range(64):
    xb = bits(x)
    for S in range(64):
        Sb = bits(S)
        parity = sum(xb[i] for i in range(6) if Sb[i]) % 2
        W_raw[x, S] = 1.0 if parity == 0 else -1.0
W = W_raw / 8.0
Q = np.round(64.0 * (W.T @ K @ W)).astype(int)

outdeg = [int(np.count_nonzero(Q[s, :])) for s in range(64)]
from collections import Counter
print("Out-degree distribution over all 64 Walsh indices:", Counter(outdeg))

indeg = [int(np.count_nonzero(Q[:, t])) for t in range(64)]
print("In-degree distribution:", Counter(indeg))

np.save('/home/claude/Q.npy', Q)
print("\nSample rows (S, out-degree, targets, weights):")
for s in [0, 1, 2, 3, 4, 5, 6, 7, 32, 43]:
    nz = np.nonzero(Q[s, :])[0]
    print(f"  S={s:2d} outdeg={outdeg[s]}: targets={nz.tolist()} weights={Q[s, nz].tolist()}")
