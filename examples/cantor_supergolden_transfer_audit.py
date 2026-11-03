"""Run the coupled Cantor--Q6--supergolden transfer audit."""
import argparse, json
from pathlib import Path
from src.cantor_supergolden_transfer import audit_cantor_q6_supergolden

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,default=Path("results/json/cantor_q6_supergolden_transfer.json"))
    args=parser.parse_args()
    report=audit_cantor_q6_supergolden()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))
    print(f"JSON report: {args.out}")

if __name__=="__main__": main()
