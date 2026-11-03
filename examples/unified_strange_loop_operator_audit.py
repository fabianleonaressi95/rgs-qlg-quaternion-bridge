import argparse,json
from pathlib import Path
from src.unified_strange_loop_operator import audit_unified_operator

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=Path("results/json/unified_strange_loop_operator.json"))
    a=p.parse_args()
    r=audit_unified_operator()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2),encoding="utf-8")
    print(json.dumps(r,indent=2))
if __name__=="__main__": main()
