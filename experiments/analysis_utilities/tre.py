# Let's verify the final validation criteria for the dataset using python code:
import pandas as pd
import numpy as np

df = pd.read_csv("data.csv")
if "norm_corrector" in df.columns:
    df["Y"] = df["norm_corrector"] / np.sqrt(df["L"])

df["u"] = df["L"] * np.sqrt(1 - df["r"])
df["A_S"] = df["Y"] * np.sqrt(df["u"])

# 1. Global exponent theta via log-log fit
log_u = np.log(df["u"])
log_Y = np.log(df["Y"])
slope, intercept = np.polyfit(log_u, log_Y, 1)
theta_global = -slope

# 2. Relative dispersion of A_S
mean_as = df["A_S"].mean()
std_as = df["A_S"].std(ddof=1)
rel_disp = (std_as / mean_as) * 100

print(f"Theta global: {theta_global:.4f} (deviation from 0.5: {abs(theta_global - 0.5):.4f})")
print(f"Mean A_S: {mean_as:.4f}, Std: {std_as:.4f}, Rel Disp: {rel_disp:.2f}%")
