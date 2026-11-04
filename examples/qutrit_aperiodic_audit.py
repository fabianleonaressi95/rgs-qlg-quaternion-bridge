"""Run the finite qutrit/Fibonacci-chain audit and print JSON."""
import json
from src.qutrit_aperiodic import audit_qutrit_aperiodic_model

if __name__ == "__main__":
    print(json.dumps(audit_qutrit_aperiodic_model(n_sites=4), indent=2))
