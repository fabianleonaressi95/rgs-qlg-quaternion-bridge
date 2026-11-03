"""Reproduce the naive/Wilson Mobius-boundary thermodynamic audit."""
import argparse, json
from pathlib import Path
from src.thermodynamics import compare_boundary_spectra

def run(input_path,output_path):
    source=Path(input_path); data=json.loads(source.read_text(encoding="utf-8"))
    report=compare_boundary_spectra(data); report["input_file"]=str(source)
    destination=Path(output_path); destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,allow_nan=False)); return report

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("input",nargs="?",default="results/json/mobius_hawking_spectrum_reference.json")
    parser.add_argument("--output",default="results/json/mobius_hawking_thermodynamics.json")
    args=parser.parse_args(); run(args.input,args.output)
