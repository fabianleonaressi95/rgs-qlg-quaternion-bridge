# Q6 exact incidence data

Source: recovered Library artifact `Q6_exact_incidence_and_H241.npz`.

SHA-256 of original NPZ: `f990699f7547fd8c26718cc2160291c4ad05c07495fe6a8f898b5e93c1265a40`.

Arrays:
- edges: (192,3), int16
- squares: (240,6), int16
- B1: (64,192), int8
- B2: (192,240), int8
- A_plaquette: (240,240), int8
- A6: (64,64), int8
- L6: (64,64), float64
- H241: (241,241), float64
- H241_eigenvalues: (241,), float64

Verified directly from the recovered NPZ:
- B1 @ B2 = zero matrix exactly
- rank(B1)=63
- rank(B2)=129
- nnz(B1)=384
- nnz(B2)=960
- entries of B1 and B2 are in {-1,0,1}.

The NPZ is the canonical source for orientation and incidence values. Do not replace it with an independently regenerated ordering. The 240 square count is combinatorial; geometric degeneracy/area is a separate property.
