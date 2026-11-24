"""Run the exact Fourier and independent matrix spectral audit for GEC-6 m=3."""
import json
from pathlib import Path

from src.gec6_a5_spectral import audit

result = audit()
out = Path("results/gec6_a5_laplacian_spectrum_m3.json")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2, sort_keys=True))
