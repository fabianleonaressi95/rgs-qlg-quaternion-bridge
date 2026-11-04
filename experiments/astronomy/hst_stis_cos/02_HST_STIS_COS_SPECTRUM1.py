#!/usr/bin/env python3

"""
CROSS-DOMAIN KPI COMPARISON ANALYZER
Unisce e confronta i KPI di Koopman estratti dai dati UV di HST (Stellari)
con i KPI di Galileo UVS (Planetari/Magnetosferici).
"""

import os
import pandas as pd
import numpy as np

# =================================================
# CONFIGURAZIONE
# =================================================

HST_KPI_PATH = "hst_koopman_kpi_results.csv"
GALILEO_KPI_PATH = "galileo_uvs_kpi_results.csv"  # Target di riferimento Galileo UVS
COMPARISON_REPORT_PATH = "cross_domain_kpi_comparison.csv"

print("="*70)
print(" CROSS-DOMAIN KPI COMPARISON ANALYZER (HST vs GALILEO UVS)")
print("="*70)

# =================================================
# 1. CARICAMENTO DEI DATASET KPI
# =================================================

if not os.path.exists(HST_KPI_PATH):
    raise RuntimeError(f"File KPI di HST non trovato: {HST_KPI_PATH}. Esegui prima l'analisi HODMD.")

df_hst = pd.read_csv(HST_KPI_PATH)
df_hst["DOMAIN"] = "HST_UV"
print(f"KPI caricati da HST: {len(df_hst)} modi analizzati.")

# Gestione opzionale/simulata di Galileo UVS se il file non è ancora presente localmente
if os.path.exists(GALILEO_KPI_PATH):
    df_galileo = pd.read_csv(GALILEO_KPI_PATH)
    df_galileo["DOMAIN"] = "GALILEO_UVS"
    print(f"KPI caricati da Galileo UVS: {len(df_galileo)} modi analizzati.")
    df_combined = pd.concat([df_hst, df_galileo], ignore_index=True)
else:
    print(f"\n[Nota] File Galileo UVS ({GALILEO_KPI_PATH}) non trovato nella working directory.")
    print("Generazione del report basato sui KPI di HST e strutturazione del framework di confronto incrociato...")
    df_combined = df_hst.copy()

# =================================================
# 2. NORMALIZZAZIONE E SINTESI STATISTICA DEI KPI
# =================================================

print("\n[2/2] Elaborazione metriche comparative di stabilità e frequenza...")

summary_stats = []

for domain, group in df_combined.groupby("DOMAIN"):
    summary_stats.append({
        "DOMAIN": domain,
        "TOTAL_MODES": len(group),
        "MEAN_MAGNITUDE": group["magnitude"].mean(),
        "STD_MAGNITUDE": group["magnitude"].std() if len(group) > 1 else 0.0,
        "MEAN_STABILITY": group["stability_metric"].mean(),
        "MAX_KPI_SCORE": group["koopman_kpi_score"].max(),
        "MEAN_KPI_SCORE": group["koopman_kpi_score"].mean()
    })

df_summary = pd.DataFrame(summary_stats)
df_summary.to_csv(COMPARISON_REPORT_PATH, index=False)

print("\n" + "="*50)
print(" REPORT COMPARATIVO CROSS-DOMAIN")
print("="*50)
print(df_summary.to_string(index=False))
print(f"\nReport globale salvato in: {COMPARISON_REPORT_PATH}")
print("="*50)
