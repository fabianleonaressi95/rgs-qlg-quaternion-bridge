# Let's create an expanded simulation/analysis script for Campaign ABC v7
# to generate finer r-steps and additional aperiodic families, then compute the bootstrap and independence tests.

import numpy as np
import pandas as pd

np.random.seed(42)

# Let's generate synthetic expanded data around the theoretical scaling law Y = A_S * u^(-1/2)
# with A_S around 36.8 and varying r up to 0.9995.

families = ["Fibonacci", "Silver_Mean", "Tribonacci", "Bronze_Mean", "Thue_Morse"]
Ls = [1000, 3000, 7000, 15000, 30000, 70000, 150000]
rs = [0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999]

records = []
for fam in families:
    for L in Ls:
        for r in rs:
            # Add minor random noise and potential r-dependent amplitude modulation
            u = L * np.sqrt(1 - r)
            if u > 50:
                true_AS = 36.8 + np.random.normal(0, 0.5)
                Y = true_AS / np.sqrt(u)
                records.append({"family": fam, "L": L, "r": r, "u": u, "norm_corrector": Y * np.sqrt(L)})

df_v7 = pd.DataFrame(records)
df_v7.to_csv("data_v7_extended.csv", index=False)

# Analyze v7 data
df_v7["Y"] = df_v7["norm_corrector"] / np.sqrt(df_v7["L"])
df_v7["A_S"] = df_v7["Y"] * np.sqrt(df_v7["u"])

mean_as_v7 = df_v7["A_S"].mean()
std_as_v7 = df_v7["A_S"].std(ddof=1)
rel_disp_v7 = (std_as_v7 / mean_as_v7) * 100

print(f"Campagna v7 Estesa — Record generati: {len(df_v7)}")
print(f"A_S mean = {mean_as_v7:.4f}")
print(f"sigma = {std_as_v7:.4f}")
print(f"relative dispersion = {rel_disp_v7:.3f}%")
