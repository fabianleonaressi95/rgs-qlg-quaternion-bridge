#!/usr/bin/env python3

import csv
from pathlib import Path


EDGE_FILE = Path(
    "rnr_hybrid_results/rnr_hybrid_edges.csv"
)

OUTPUT_DIR = Path(
    "rnr_experimental_calibration"
)


# Dati ricavati dalla letteratura.
#
# Non sono cinque misure indipendenti:
# alcuni sono intervalli di trasferimento,
# altri sono ordini di grandezza o limiti inferiori.
EXPERIMENTAL_DATA = {
    ("Y122", "W48"): {
        "rate": "",
        "uncertainty": "",
        "kie": "",
        "temperature": "",
        "status": "not separately resolved",
        "note": (
            "No edge-specific experimental rate "
            "isolated for Y122-W48."
        ),
    },
    ("W48", "Y356"): {
        "rate": "",
        "uncertainty": "",
        "kie": "",
        "temperature": "",
        "status": "not separately resolved",
        "note": (
            "No edge-specific experimental rate "
            "isolated for W48-Y356."
        ),
    },
    ("Y356", "Y731"): {
        "rate": "6.0e3-1.6e4",
        "uncertainty": "",
        "kie": "",
        "temperature": "",
        "status": "interfacial RT range",
        "note": (
            "Conformationally dependent radical "
            "transfer across the alpha/beta interface."
        ),
    },
    ("Y731", "Y730"): {
        "rate": "1.0e5",
        "uncertainty": "",
        "kie": "",
        "temperature": "",
        "status": "PCET step order",
        "note": (
            "Order of magnitude reported for "
            "photoinitiated PCET studies."
        ),
    },
    ("Y730", "C439"): {
        "rate": "1.4e4",
        "uncertainty": "0.4e4",
        "kie": "",
        "temperature": "",
        "status": "lower bound",
        "note": (
            "Lower bound for the substrate C-H-limited "
            "forward PCET step; not an isolated "
            "Y730-C439 rate."
        ),
    },
}


def read_edges():
    if not EDGE_FILE.exists():
        raise FileNotFoundError(
            f"File non trovato: {EDGE_FILE}\n"
            "Esegui lo script dalla cartella principale "
            "del progetto."
        )

    with EDGE_FILE.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_calibrated_edges(rows):
    OUTPUT_DIR.mkdir(exist_ok=True)

    output = (
        OUTPUT_DIR
        / "rnr_experimental_calibration.csv"
    )

    columns = [
        "source",
        "target",
        "source_chain",
        "target_chain",
        "geometry_source",
        "d_ET_geometry_A",
        "k_exp_s-1_or_range",
        "uncertainty_s-1",
        "KIE_exp",
        "temperature_K",
        "evidence_status",
        "source_note",
    ]

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)

        for row in rows:
            source = row.get("source", "")
            target = row.get("target", "")

            data = EXPERIMENTAL_DATA.get(
                (source, target),
                {
                    "rate": "",
                    "uncertainty": "",
                    "kie": "",
                    "temperature": "",
                    "status": "no datum",
                    "note": "",
                },
            )

            writer.writerow([
                source,
                target,
                row.get("source_chain", ""),
                row.get("target_chain", ""),
                row.get("geometry_source", ""),
                row.get("centroid_distance_A", ""),
                data["rate"],
                data["uncertainty"],
                data["kie"],
                data["temperature"],
                data["status"],
                data["note"],
            ])

    return output


def write_global_controls():
    output = (
        OUTPUT_DIR
        / "rnr_global_experimental_controls.csv"
    )

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)

        writer.writerow([
            "quantity",
            "value",
            "unit",
            "temperature_K",
            "scope",
            "status",
        ])

        writer.writerow([
            "KIE",
            "",
            "dimensionless",
            "",
            "global/process-level",
            (
                "No single edge-resolved H/D value "
                "inserted automatically."
            ),
        ])

        writer.writerow([
            "k_C-H_lower_bound",
            "1.4e4",
            "s^-1",
            "",
            "substrate C-H-limited step",
            (
                "Lower bound; not assigned to every "
                "pathway edge."
            ),
        ])

        writer.writerow([
            "interfacial_radical_transfer",
            "0.6e4-1.6e4",
            "s^-1",
            "",
            "Y356-Y731 interface region",
            (
                "Conformationally dependent range."
            ),
        ])

        writer.writerow([
            "PCET_step_order",
            "1.0e5",
            "s^-1",
            "",
            "photoinitiated PCET",
            (
                "Order of magnitude; not a direct "
                "measurement of one isolated edge."
            ),
        ])

    return output


def print_summary(rows):
    print("Calibrazione sperimentale RNR")
    print()
    print(f"Archi letti: {len(rows)}")
    print()

    for row in rows:
        source = row.get("source", "")
        target = row.get("target", "")
        data = EXPERIMENTAL_DATA.get(
            (source, target),
            {},
        )

        print(
            f"{source} -> {target}: "
            f"{data.get('status', 'no datum')}"
        )

        if data.get("rate"):
            print(
                f"  k sperimentale: "
                f"{data['rate']} s^-1"
            )
        else:
            print(
                "  k sperimentale: "
                "non risolto per questo arco"
            )


def main():
    rows = read_edges()

    calibrated_file = write_calibrated_edges(rows)
    controls_file = write_global_controls()

    print_summary(rows)

    print()
    print("File prodotti:")
    print(f"  {calibrated_file}")
    print(f"  {controls_file}")


if __name__ == "__main__":
    main()
