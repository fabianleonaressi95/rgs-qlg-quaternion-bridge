import pandas as pd
import numpy as np
import statsmodels.api as sm

# Load data
df = pd.read_csv("data_v7_extended.csv")

# Ensure Y and u are correctly computed
if "Y" not in df.columns and "norm_corrector" in df.columns:
    df["Y"] = df["norm_corrector"] / np.sqrt(df["L"])

df["u"] = df["L"] * np.sqrt(1 - df["r"])

# Asymptotic filter
df_filtered = df[df["u"] > 100].copy()

print("N punti filtrati:", len(df_filtered))

# Separable model fit
X1 = np.log(df_filtered["L"])
X2 = np.log(1 - df_filtered["r"])
Ylog = np.log(df_filtered["Y"])

X = sm.add_constant(np.column_stack([X1, X2]))
model = sm.OLS(Ylog, X).fit()

a, alpha_neg, beta_neg = model.params
alpha = -alpha_neg
beta = -beta_neg

print("\n--- SEPARABLE SCALING FIT ---")
print(f"A = {np.exp(a):.6f}")
print(f"alpha = {alpha:.6f}")
print(f"beta = {beta:.6f}")
print(f"R2 sep = {model.rsquared:.6f}")

# u-scaling model fit
Xu = sm.add_constant(np.log(df_filtered["u"]))
model_u = sm.OLS(Ylog, Xu).fit()
A_u = model_u.params[0]
theta = -model_u.params[1]

print("\n--- u SCALING MODEL ---")
print(f"A = {np.exp(A_u):.6f}")
print(f"theta = {theta:.6f}")
print(f"R2 u = {model_u.rsquared:.6f}")

# Bootstrap
Nboot = 2000
results = []
for i in range(Nboot):
    sample = df_filtered.sample(len(df_filtered), replace=True)
    Xb = sm.add_constant(np.column_stack([np.log(sample["L"]), np.log(1 - sample["r"])]))
    yb = np.log(sample["Y"])
    fit = sm.OLS(yb, Xb).fit()
    results.append([-fit.params[1], -fit.params[2]])

results = np.array(results)
alpha_ci = np.percentile(results[:, 0], [2.5, 97.5])
beta_ci = np.percentile(results[:, 1], [2.5, 97.5])

print("\n--- BOOTSTRAP 95% ---")
print(f"alpha CI: {alpha_ci}")
print(f"beta CI: {beta_ci}")

# Amplitude
df_filtered["A_S"] = df_filtered["Y"] * np.sqrt(df_filtered["u"])
mean_as = df_filtered["A_S"].mean()
std_as = df_filtered["A_S"].std()
disp = 100 * std_as / mean_as

print("\n--- AMPLITUDE ---")
print(f"Mean A_S: {mean_as:.4f}")
print(f"Sigma: {std_as:.4f}")
print(f"Dispersion: {disp:.2f}%")

df_filtered.to_csv("results_campaign_v8.csv", index=False)
