from src.q6_prism_symmetry_obstruction import audit
import json

if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
