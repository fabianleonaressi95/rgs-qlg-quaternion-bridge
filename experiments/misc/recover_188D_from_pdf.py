# -*- coding: utf-8 -*-
# 188D / 64B — forensic recovery from the audit PDF.
# Does NOT invent X_Oh. It extracts the embedded notebook source,
# searches for the original geometry construction, and dumps exact
# source contexts for reconstruction.

import os, re, json, hashlib
from pathlib import Path

OUT = Path("/kaggle/working/188D_RECOVERY")
OUT.mkdir(parents=True, exist_ok=True)

# 1) Find audit PDF
pdfs = []
for base in (Path("/kaggle/input"), Path("/kaggle/working"), Path("/kaggle/temp")):
    if base.exists():
        pdfs += list(base.rglob("*.pdf"))

preferred = [x for x in pdfs if x.name.lower() in ("formule.pdf", "formule(1).pdf")]
if preferred:
    pdf_path = preferred[0]
elif pdfs:
    pdf_path = max(pdfs, key=lambda x: x.stat().st_size)
else:
    raise RuntimeError("Nessun PDF trovato.")

print("PDF:", pdf_path)

# 2) Extract text
try:
    import fitz
    doc = fitz.open(str(pdf_path))
    text = "\n".join(page.get_text("text") for page in doc)
    backend = "PyMuPDF"
except Exception:
    from pypdf import PdfReader
    reader = PdfReader(str(pdf_path))
    text = "\n".join((page.extract_text() or "") for page in reader.pages)
    backend = "pypdf"

(OUT / "pdf_text_raw.txt").write_text(text, encoding="utf-8")

# 3) Locate embedded notebook source
m = re.search(
    r"--- FILE:\s*(?:/kaggle/working/)?\.virtual_documents/__notebook_source__\.ipynb\s*---",
    text, re.I
)
if not m:
    raise RuntimeError(
        "Il PDF non contiene il blocco __notebook_source__.ipynb riconoscibile. "
        "Il testo grezzo e' comunque in pdf_text_raw.txt."
    )

rest = text[m.end():]
m2 = re.search(r"\n--- FILE:", rest, re.I)
raw_source = rest[:m2.start()] if m2 else rest

# 4) Remove PDF line numbers such as '206: code'
lines = []
for line in raw_source.splitlines():
    q = re.match(r"^\s*\d+:\s?(.*)$", line)
    lines.append(q.group(1) if q else line)

source = "\n".join(lines)
(OUT / "notebook_source_recovered.py").write_text(source, encoding="utf-8")

# 5) Forensic search
patterns = [
    r"\bX_Oh\s*=",
    r"\bX\s*=.*64",
    r"64B[_ ]?Oh",
    r"reconstruct.*X",
    r"construct.*X",
    r"build.*X",
    r"generate.*X",
    r"sqrt\(5\)",
    r"phi",
    r"vertices",
    r"coordinates",
]

hits = {}
contexts = []

for pat in patterns:
    ms = list(re.finditer(pat, source, re.I))
    hits[pat] = len(ms)
    for mm in ms[:20]:
        a = max(0, mm.start()-1800)
        b = min(len(source), mm.end()+3500)
        contexts.append({
            "pattern": pat,
            "offset": mm.start(),
            "context": source[a:b]
        })

(OUT / "X_Oh_contexts.json").write_text(
    json.dumps(contexts, indent=2, ensure_ascii=False),
    encoding="utf-8"
)

# 6) Core pipeline targets
targets = [
    "def coulomb_hessian",
    "def tangent_basis",
    "def constrained_quantities",
    "lambda_force",
    "def build_field_operator",
    "K_field",
    "B_T2U_TAN",
    "J_A",
    "J_B",
    "P_G",
    "Vmass",
    "logm",
    "Omega",
    "D = -Omega",
    "Kinv",
]

target_hits = {}
snippets = []

for target in targets:
    ms = list(re.finditer(re.escape(target), source, re.I))
    target_hits[target] = len(ms)
    if ms:
        mm = ms[0]
        a = max(0, mm.start()-500)
        b = min(len(source), mm.start()+10000)
        snippets.append(
            "\n" + "="*100 +
            "\nTARGET: " + target +
            "\n" + "="*100 +
            "\n" + source[a:b]
        )

(OUT / "core_source_snippets.txt").write_text(
    "\n".join(snippets), encoding="utf-8"
)

# 7) References to geometry files in the PDF
refs = sorted(set(re.findall(
    r"[^\s\"']*(?:X_Oh|64B_Oh|reconstructed)[^\s\"']*\.npy",
    text, re.I
)))

manifest = {
    "pdf": str(pdf_path),
    "pdf_sha256": hashlib.sha256(pdf_path.read_bytes()).hexdigest(),
    "backend": backend,
    "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
    "source_chars": len(source),
    "xoh_contexts": len(contexts),
    "geometry_file_refs": refs,
    "target_hits": target_hits,
}

(OUT / "recovery_manifest.json").write_text(
    json.dumps(manifest, indent=2, ensure_ascii=False),
    encoding="utf-8"
)

print("\n" + "="*80)
print("188D FORENSIC RECOVERY")
print("="*80)
print("Output:", OUT)
print("Notebook source:", OUT / "notebook_source_recovered.py")
print("X_Oh contexts:", len(contexts))
print("\nTarget hits:")
for k, v in target_hits.items():
    print(f"  {k:28s} {v}")

print("\nGeometry references found in PDF:")
for x in refs:
    print(" ", x)

print("\nNON INVENTO X_Oh.")
print("Il prossimo passaggio puo' usare questi artefatti per estrarre")
print("la costruzione originale e ricostruire la pipeline completa.")
