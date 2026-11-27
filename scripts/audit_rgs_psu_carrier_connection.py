#!/usr/bin/env python3
"""Run the RGS–PSU / carrier connection comparison and print its JSON report."""
import json
from src.rgs_psu_carrier_connection_audit import audit_rgs_psu_carrier_connection

if __name__ == "__main__":
    print(json.dumps(audit_rgs_psu_carrier_connection(), indent=2, sort_keys=True))
