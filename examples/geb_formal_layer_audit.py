"""Reproducible Godel--Escher--Bach audit."""
import argparse, json
from pathlib import Path
from src.geb_formal_layer import audit_geb

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=Path("results/json/geb_formal_layer.json"))
    a=p.parse_args()
    report=audit_geb()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))
    print(f"JSON report: {a.out}")

if __name__=="__main__": main()
