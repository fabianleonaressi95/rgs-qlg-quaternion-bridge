"""Run the exact finite audit for the GEC-6 trigram decomposition."""
from src.gec6_trigram_decomposition import audit
import json

if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
