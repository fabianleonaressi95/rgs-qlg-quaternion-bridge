#!/usr/bin/env python3
"""
QVSM-RNR source collector.

Obiettivi:
1. Registra fonti primarie RNR/PCET/NEO-QM-MM.
2. Scarica legalmente contenuti pubblici accessibili via URL.
3. Cerca sul DOI ACS link a Supporting Information o file dati.
4. Genera:
   - sources.csv
   - claims.csv
   - neo_raw_samples_template.csv
   - author_data_request.txt
   - manifest.json
5. Non inventa dati configurazione-per-configurazione.
"""
import numpy as np
import argparse
import csv
import hashlib
import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import urljoin, urlparse

import pandas as pd
import requests
from bs4 import BeautifulSoup

# Shared atomic writer; supports both package and direct-script execution.
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared_io import atomic_write_bytes


USER_AGENT = (
    "QVSM-RNR-reproducibility-audit/1.0 "
    "(academic provenance collection; contact: local-user)"
)

TIMEOUT_SECONDS = 45
SLEEP_SECONDS = 1.0

PRIMARY_SOURCES = [
    {
        "source_id": "neo_nqe_2024",
        "title": (
            "Nuclear Quantum Effects in Quantum Mechanical/Molecular "
            "Mechanical Free Energy Simulations of Ribonucleotide Reductase"
        ),
        "doi": "10.1021/jacs.4c13955",
        "year": 2024,
        "journal": "Journal of the American Chemical Society",
        "landing_url": "https://pubs.acs.org/doi/10.1021/jacs.4c13955",
        "pmc_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11625381/",
        "role": "NEO_QM_MM_NQE",
        "relevance": (
            "Y731-Y730 NEO-DFT/QM-MM reweighting; "
            "NQE free-energy barrier shifts; H/D framework"
        ),
        "expected_public_content": (
            "Main article; Supporting Information PDF with computational "
            "details, string convergence, bootstrapping, protonated "
            "glutamate results, basis-set benchmarking"
        ),
    },
    {
        "source_id": "nonadiabaticity_2024",
        "title": (
            "Probing Nonadiabaticity of Proton-Coupled Electron Transfer "
            "in Ribonucleotide Reductase"
        ),
        "doi": "10.1021/acs.jpclett.3c03552",
        "year": 2024,
        "journal": "Journal of Physical Chemistry Letters",
        "landing_url": "https://pubs.acs.org/doi/10.1021/acs.jpclett.3c03552",
        "pmc_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11586673/",
        "role": "PCET_REGIME_CLASSIFICATION",
        "relevance": (
            "Y356-Y731 nonadiabatic; Y731-Y730 adiabatic; "
            "Y730-C439 intermediate"
        ),
        "expected_public_content": (
            "Main article; methodology and PCET-regime classification"
        ),
    },
    {
        "source_id": "interfacial_direct_pcet_2023",
        "title": (
            "Direct Proton-Coupled Electron Transfer between Interfacial "
            "Tyrosines in Ribonucleotide Reductase"
        ),
        "doi": "10.1021/jacs.2c13615",
        "year": 2023,
        "journal": "Journal of the American Chemical Society",
        "landing_url": "https://pubs.acs.org/doi/10.1021/jacs.2c13615",
        "pmc_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10344599/",
        "role": "INTERFACIAL_PCET_QM_MM",
        "relevance": (
            "Y356-Y731 direct PCET; comparison with water-mediated "
            "mechanism; role of flipped-out conformation"
        ),
        "expected_public_content": (
            "Main article; QM/MM pathway comparison"
        ),
    },
    {
        "source_id": "water_network_2020",
        "title": (
            "Conformational Motions and Water Networks at the α/β "
            "Interface in E. coli Ribonucleotide Reductase"
        ),
        "doi": "10.1021/jacs.0c04325",
        "year": 2020,
        "journal": "Journal of the American Chemical Society",
        "landing_url": "https://pubs.acs.org/doi/10.1021/jacs.0c04325",
        "pmc_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7594210/",
        "role": "WATER_NETWORK_MD",
        "relevance": (
            "Y731 flipped-out and stacked states; interfacial water; "
            "bridging-water configurations"
        ),
        "expected_public_content": (
            "Main article; MD/free-energy simulation analysis"
        ),
    },
    {
        "source_id": "oxygen_endor_2021",
        "title": (
            "Detection of Water Molecules on the Radical Transfer Pathway "
            "of Ribonucleotide Reductase by 17O Electron-Nuclear Double "
            "Resonance Spectroscopy"
        ),
        "doi": "10.1021/jacs.1c02020",
        "year": 2021,
        "journal": "Journal of the American Chemical Society",
        "landing_url": "https://pubs.acs.org/doi/10.1021/jacs.1c02020",
        "pmc_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8154519/",
        "role": "WATER_ENDOR_EXPERIMENT",
        "relevance": (
            "17O-ENDOR evidence for water associated with radical "
            "transfer pathway intermediates"
        ),
        "expected_public_content": (
            "Main article; spectroscopic evidence and analysis"
        ),
    },
    {
        "source_id": "hbond_endor_2015",
        "title": (
            "Hydrogen Bond Network between Amino Acid Radical "
            "Intermediates on the Proton-Coupled Electron Transfer "
            "Pathway of E. coli α2 Ribonucleotide Reductase"
        ),
        "doi": "10.1021/ja510513z",
        "year": 2015,
        "journal": "Journal of the American Chemical Society",
        "landing_url": "https://pubs.acs.org/doi/10.1021/ja510513z",
        "pmc_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4516058/",
        "role": "HBOND_ENDOR_EXPERIMENT",
        "relevance": (
            "ENDOR/DFT evidence for H-bond network around Y731/Y730/C439"
        ),
        "expected_public_content": (
            "Main article; ENDOR and DFT assignments"
        ),
    },
    {
        "source_id": "structure_6w4x_2020",
        "title": (
            "Structure of a Trapped Radical Transfer Pathway within a "
            "Ribonucleotide Reductase Holocomplex"
        ),
        "doi": "10.1126/science.aba6794",
        "year": 2020,
        "journal": "Science",
        "landing_url": "https://www.science.org/doi/10.1126/science.aba6794",
        "pmc_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7774503/",
        "role": "STRUCTURE_6W4X",
        "relevance": (
            "Cryo-EM holocomplex; structural basis of the radical transfer "
            "pathway; PDB 6W4X"
        ),
        "expected_public_content": (
            "Main article; structural model; PDB provenance"
        ),
    },
]


CLAIMS = [
    {
        "claim_id": "neo_method",
        "source_id": "neo_nqe_2024",
        "claim": (
            "NEO-DFT treats selected protonic nuclei quantum mechanically "
            "with electrons and is combined with QM/MM finite-temperature "
            "string/umbrella sampling via reweighting."
        ),
        "evidence_type": "QM_DERIVED_METHOD",
        "scope": "METHOD",
        "status": "SUPPORTED",
    },
    {
        "claim_id": "neo_nqe_sign",
        "source_id": "neo_nqe_2024",
        "claim": (
            "Nuclear quantum effects can either raise or lower the free "
            "energy barrier in the studied RNR PCET scenarios."
        ),
        "evidence_type": "QM_DERIVED_RESULT",
        "scope": "Y731_Y730_PATHWAY_SCENARIOS",
        "status": "SUPPORTED",
    },
    {
        "claim_id": "e4_regime",
        "source_id": "nonadiabaticity_2024",
        "claim": (
            "Y731-Y730 is classified as vibronically adiabatic, whereas "
            "Y356-Y731 is nonadiabatic and Y730-C439 intermediate."
        ),
        "evidence_type": "QM_DERIVED_REGIME",
        "scope": "EDGE_REGIME",
        "status": "SUPPORTED",
    },
    {
        "claim_id": "y731_conformations",
        "source_id": "water_network_2020",
        "claim": (
            "Y731 samples flipped-out and stacked conformations that "
            "support different portions of the PCET pathway."
        ),
        "evidence_type": "MD_FREE_ENERGY",
        "scope": "Y731_CONFORMATION",
        "status": "SUPPORTED",
    },
    {
        "claim_id": "interfacial_water",
        "source_id": "water_network_2020",
        "claim": (
            "Interfacial water and hydrogen-bonded water chains modulate "
            "the local PCET environment."
        ),
        "evidence_type": "MD_QM_MM_SUPPORTED",
        "scope": "INTERFACE",
        "status": "SUPPORTED",
    },
    {
        "claim_id": "water_endor",
        "source_id": "oxygen_endor_2021",
        "claim": (
            "17O-ENDOR provides spectroscopic evidence for water molecules "
            "associated with radical-transfer pathway intermediates."
        ),
        "evidence_type": "MEASURED_SPECTROSCOPY",
        "scope": "WATER_NETWORK",
        "status": "SUPPORTED",
    },
    {
        "claim_id": "qed_water",
        "source_id": "neo_nqe_2024",
        "claim": (
            "Macroscopic QED water coherence domains drive RNR PCET."
        ),
        "evidence_type": "NONE",
        "scope": "QED_WATER",
        "status": "NOT_SUPPORTED_BY_RNR_SOURCE",
    },
]


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def safe_filename(text, max_length=120):
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    return text[:max_length].strip("_")


def get_session():
    session = requests.Session()
    session.headers.update({
        "User-Agent": USER_AGENT,
        "Accept": (
            "text/html,application/pdf,application/xml,"
            "application/json;q=0.9,*/*;q=0.8"
        ),
    })
    return session


def request_url(session, url):
    try:
        response = session.get(
            url,
            timeout=TIMEOUT_SECONDS,
            allow_redirects=True,
        )
        return response
    except requests.RequestException as error:
        return None, str(error)


def save_response_content(response, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_bytes(response.content, destination)

    return {
        "path": str(destination),
        "bytes": len(response.content),
        "sha256": sha256_file(destination),
        "content_type": response.headers.get("Content-Type", ""),
        "final_url": response.url,
        "status_code": response.status_code,
    }


def is_candidate_download_url(url):
    lowered = url.lower()

    markers = [
        "suppl",
        "support",
        "si_",
        "media",
        "download",
        ".pdf",
        ".zip",
        ".csv",
        ".xlsx",
        ".xls",
        ".json",
        ".txt",
        ".dat",
    ]

    return any(marker in lowered for marker in markers)


def extract_links(base_url, html):
    soup = BeautifulSoup(html, "html.parser")
    links = []

    for anchor in soup.find_all("a", href=True):
        href = anchor.get("href", "").strip()
        text = " ".join(anchor.get_text(" ", strip=True).split())

        if not href:
            continue

        absolute = urljoin(base_url, href)

        links.append({
            "text": text,
            "url": absolute,
            "candidate_download": is_candidate_download_url(absolute)
                                  or is_candidate_download_url(text),
        })

    unique = []
    seen = set()

    for link in links:
        key = (link["text"], link["url"])

        if key not in seen:
            unique.append(link)
            seen.add(key)

    return unique


def download_public_url(session, url, destination):
    result = {
        "url": url,
        "download_status": "NOT_ATTEMPTED",
        "error": "",
        "path": "",
        "bytes": np.nan,
        "sha256": "",
        "content_type": "",
        "final_url": "",
        "status_code": np.nan,
    }

    try:
        response = session.get(
            url,
            timeout=TIMEOUT_SECONDS,
            allow_redirects=True,
        )

        result["final_url"] = response.url
        result["status_code"] = response.status_code
        result["content_type"] = response.headers.get("Content-Type", "")

        if response.status_code != 200:
            result["download_status"] = "HTTP_ERROR"
            result["error"] = f"HTTP {response.status_code}"
            return result

        content_type = response.headers.get("Content-Type", "").lower()

        if (
            "text/html" in content_type
            and destination.suffix.lower() not in {".html", ".htm"}
        ):
            result["download_status"] = "HTML_NOT_FILE"
            result["error"] = (
                "URL returned HTML rather than a directly downloadable file."
            )
            return result

        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(response.content)

        result["download_status"] = "DOWNLOADED"
        result["path"] = str(destination)
        result["bytes"] = len(response.content)
        result["sha256"] = sha256_file(destination)

        return result

    except requests.RequestException as error:
        result["download_status"] = "REQUEST_ERROR"
        result["error"] = str(error)
        return result


def source_record_to_csv_row(source):
    return {
        "source_id": source["source_id"],
        "title": source["title"],
        "doi": source["doi"],
        "year": source["year"],
        "journal": source["journal"],
        "role": source["role"],
        "relevance": source["relevance"],
        "landing_url": source["landing_url"],
        "pmc_url": source["pmc_url"],
        "expected_public_content": source["expected_public_content"],
    }


def make_raw_template(path):
    columns = [
        "sample_id",
        "system",
        "edge",
        "configuration_id",
        "window_id",
        "string_image_id",
        "R_OO_A",
        "delta_A",
        "weight_classical",
        "E_DFT_kcal_mol",
        "E_NEO_H_kcal_mol",
        "E_NEO_D_kcal_mol",
        "temperature_K",
        "functional",
        "electronic_basis",
        "protonic_basis",
        "QM_region",
        "MM_environment",
        "source_reference",
        "source_dataset",
        "evidence_status",
        "notes",
    ]

    example = {
        "sample_id": "REPLACE_WITH_REAL_SAMPLE_ID",
        "system": "alpha_E623_mediated_or_alpha_prime_direct",
        "edge": "e4",
        "configuration_id": "REPLACE_WITH_FRAME_OR_SNAPSHOT_ID",
        "window_id": "REPLACE_WITH_UMBRELLA_WINDOW",
        "string_image_id": "REPLACE_WITH_STRING_IMAGE",
        "R_OO_A": "NA",
        "delta_A": "NA",
        "weight_classical": "NA",
        "E_DFT_kcal_mol": "NA",
        "E_NEO_H_kcal_mol": "NA",
        "E_NEO_D_kcal_mol": "NA",
        "temperature_K": "298.15",
        "functional": "REPLACE_WITH_PUBLISHED_METHOD",
        "electronic_basis": "REPLACE_WITH_PUBLISHED_METHOD",
        "protonic_basis": "REPLACE_WITH_PUBLISHED_METHOD",
        "QM_region": "REPLACE_WITH_QM_REGION",
        "MM_environment": "REPLACE_WITH_MM_SETUP",
        "source_reference": "10.1021/jacs.4c13955",
        "source_dataset": "SI_OR_AUTHOR_SUPPLIED_ARCHIVE",
        "evidence_status": "QM_DERIVED",
        "notes": (
            "Do not populate with invented values. One row per actual "
            "QM/MM configuration evaluated by both DFT and NEO-DFT."
        ),
    }

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerow(example)


def make_surface_template(path):
    columns = [
        "edge",
        "system",
        "isotope",
        "reactant_basin_id",
        "transition_state_id",
        "deltaG_classical_kcal_mol",
        "deltaG_quantum_kcal_mol",
        "Q_deltaDeltaG_NQE_kcal_mol",
        "temperature_K",
        "method",
        "source_reference",
        "evidence_scope",
        "evidence_status",
        "notes",
    ]

    rows = [
        {
            "edge": "e4",
            "system": "alpha_E623_mediated",
            "isotope": "H",
            "reactant_basin_id": "published_pathway",
            "transition_state_id": "published_pathway",
            "deltaG_classical_kcal_mol": "4.2",
            "deltaG_quantum_kcal_mol": "5.0",
            "Q_deltaDeltaG_NQE_kcal_mol": "0.8",
            "temperature_K": "298.15",
            "method": "NEO-DFT/QM-MM reweighting",
            "source_reference": "10.1021/jacs.4c13955",
            "evidence_scope": "MODEL_PATHWAY_SCENARIO",
            "evidence_status": "PUBLISHED_QM_DERIVED",
            "notes": (
                "Published H-pathway scenario. Small shift subject to "
                "chemical-accuracy caution. Not experimental edge rate."
            ),
        },
        {
            "edge": "e4",
            "system": "alpha_prime_direct",
            "isotope": "H",
            "reactant_basin_id": "published_pathway",
            "transition_state_id": "published_pathway",
            "deltaG_classical_kcal_mol": "4.2",
            "deltaG_quantum_kcal_mol": "-0.2",
            "Q_deltaDeltaG_NQE_kcal_mol": "-4.4",
            "temperature_K": "298.15",
            "method": "NEO-DFT/QM-MM reweighting",
            "source_reference": "10.1021/jacs.4c13955",
            "evidence_scope": "MODEL_PATHWAY_SCENARIO",
            "evidence_status": "PUBLISHED_QM_DERIVED",
            "notes": (
                "Published direct alpha-prime pathway scenario. Not to be "
                "merged with pre-turnover alpha scenario."
            ),
        },
        {
            "edge": "e4",
            "system": "alpha_E623_mediated",
            "isotope": "D",
            "reactant_basin_id": "REQUIRES_ISOTOPE_RESOLVED_DATA",
            "transition_state_id": "REQUIRES_ISOTOPE_RESOLVED_DATA",
            "deltaG_classical_kcal_mol": "NA",
            "deltaG_quantum_kcal_mol": "NA",
            "Q_deltaDeltaG_NQE_kcal_mol": "NA",
            "temperature_K": "298.15",
            "method": "NEO-DFT/QM-MM reweighting",
            "source_reference": "10.1021/jacs.4c13955",
            "evidence_scope": "UNAVAILABLE",
            "evidence_status": "NEEDS_RAW_OR_SI_EXTRACTION",
            "notes": (
                "Do not infer D barrier or KIE from H-only published shift."
            ),
        },
    ]

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def make_author_request(path):
    text = """Subject: Request for configuration-level NEO-QM/MM reweighting data for E. coli RNR

Dear Dr. Chow, Dr. Reinhardt, and Dr. Hammes-Schiffer,

I am building a reproducible audit-oriented coarse-grained representation
of the Y731–Y730 proton-coupled electron transfer step in E. coli
ribonucleotide reductase, using your study:

Nuclear Quantum Effects in Quantum Mechanical/Molecular Mechanical Free
Energy Simulations of Ribonucleotide Reductase
DOI: 10.1021/jacs.4c13955

I would like to reproduce the NEO-DFT reweighting at the level of the
published free-energy surfaces, without changing the physical model.

If available, could you please indicate a public archive or share the
configuration-level data required for the reweighting? The most useful
fields would be:

1. sampled QM/MM configuration identifiers and/or coordinates;
2. alpha versus alpha-prime ensemble assignment;
3. finite-temperature string/umbrella window or string-image identity;
4. collective variables R_OO and proton-transfer coordinate delta;
5. classical umbrella/WHAM/MBAR weights or equivalent reweighting weights;
6. E_DFT for each evaluated configuration;
7. E_NEO,H and E_NEO,D, or equivalently deltaE_NQE,H and deltaE_NQE,D;
8. temperature and definitions of reactant and transition-state basins;
9. exact QM region, electronic functional/basis, protonic basis, and MM setup;
10. numerical free-energy grids if the configuration-level data cannot be shared.

I will retain the provenance, cite the article and any associated dataset,
and clearly label these values as QM-derived rather than experimental.

Best regards,
[YOUR NAME]
"""
    path.write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Collect public RNR/NEO-QM-MM sources and create "
            "reproducible data templates."
        )
    )

    parser.add_argument(
        "--outdir",
        default="rnr_neo_sources",
        help="Directory output.",
    )

    parser.add_argument(
        "--download-pmc-html",
        action="store_true",
        help="Download public PMC HTML pages.",
    )

    parser.add_argument(
        "--download-acs-candidates",
        action="store_true",
        help=(
            "Attempt download of publicly reachable ACS candidate links "
            "identified from DOI landing pages."
        ),
    )

    parser.add_argument(
        "--max-acs-candidates-per-source",
        type=int,
        default=20,
        help="Max candidate download links checked per ACS source.",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    html_dir = outdir / "html"
    downloads_dir = outdir / "downloads"
    tables_dir = outdir / "tables"
    templates_dir = outdir / "templates"
    logs_dir = outdir / "logs"

    for directory in [
        html_dir,
        downloads_dir,
        tables_dir,
        templates_dir,
        logs_dir,
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    session = get_session()

    source_rows = []
    link_rows = []
    download_rows = []

    for source in PRIMARY_SOURCES:
        print(f"\nSOURCE: {source['source_id']}")
        source_rows.append(source_record_to_csv_row(source))

        if args.download_pmc_html:
            try:
                response = session.get(
                    source["pmc_url"],
                    timeout=TIMEOUT_SECONDS,
                    allow_redirects=True,
                )

                pmc_file = (
                    html_dir
                    / f"{source['source_id']}_PMC.html"
                )

                if response.status_code == 200:
                    atomic_write_bytes(response.content, pmc_file)

                    download_rows.append({
                        "source_id": source["source_id"],
                        "kind": "PMC_HTML",
                        "original_url": source["pmc_url"],
                        "final_url": response.url,
                        "status": "DOWNLOADED",
                        "http_status": response.status_code,
                        "content_type": response.headers.get(
                            "Content-Type",
                            "",
                        ),
                        "path": str(pmc_file),
                        "bytes": len(response.content),
                        "sha256": sha256_file(pmc_file),
                        "error": "",
                    })

                    print(f"  PMC HTML downloaded: {pmc_file.name}")

                else:
                    download_rows.append({
                        "source_id": source["source_id"],
                        "kind": "PMC_HTML",
                        "original_url": source["pmc_url"],
                        "final_url": response.url,
                        "status": "HTTP_ERROR",
                        "http_status": response.status_code,
                        "content_type": response.headers.get(
                            "Content-Type",
                            "",
                        ),
                        "path": "",
                        "bytes": np.nan,
                        "sha256": "",
                        "error": f"HTTP {response.status_code}",
                    })

            except requests.RequestException as error:
                download_rows.append({
                    "source_id": source["source_id"],
                    "kind": "PMC_HTML",
                    "original_url": source["pmc_url"],
                    "final_url": "",
                    "status": "REQUEST_ERROR",
                    "http_status": np.nan,
                    "content_type": "",
                    "path": "",
                    "bytes": np.nan,
                    "sha256": "",
                    "error": str(error),
                })

        # Scarica HTML landing ACS per trovare SI o file dataset.
        try:
            response = session.get(
                source["landing_url"],
                timeout=TIMEOUT_SECONDS,
                allow_redirects=True,
            )

            acs_html_file = (
                html_dir
                / f"{source['source_id']}_landing.html"
            )

            if response.status_code == 200:
                atomic_write_bytes(response.content, acs_html_file)

                download_rows.append({
                    "source_id": source["source_id"],
                    "kind": "LANDING_HTML",
                    "original_url": source["landing_url"],
                    "final_url": response.url,
                    "status": "DOWNLOADED",
                    "http_status": response.status_code,
                    "content_type": response.headers.get(
                        "Content-Type",
                        "",
                    ),
                    "path": str(acs_html_file),
                    "bytes": len(response.content),
                    "sha256": sha256_file(acs_html_file),
                    "error": "",
                })

                links = extract_links(response.url, response.text)

                candidates = [
                    link for link in links
                    if link["candidate_download"]
                ]

                for link in candidates:
                    link_rows.append({
                        "source_id": source["source_id"],
                        "link_text": link["text"],
                        "url": link["url"],
                        "candidate_download": True,
                    })

                print(
                    f"  landing downloaded; "
                    f"candidate links found: {len(candidates)}"
                )

                if args.download_acs_candidates:
                    candidates = candidates[
                        :args.max_acs_candidates_per_source
                    ]

                    for index, link in enumerate(candidates, start=1):
                        suffix = Path(
                            urlparse(link["url"]).path
                        ).suffix.lower()

                        if suffix not in {
                            ".pdf",
                            ".zip",
                            ".csv",
                            ".xlsx",
                            ".xls",
                            ".json",
                            ".txt",
                            ".dat",
                        }:
                            suffix = ".bin"

                        filename = (
                            f"{source['source_id']}_candidate_"
                            f"{index:02d}{suffix}"
                        )

                        destination = downloads_dir / filename

                        result = download_public_url(
                            session,
                            link["url"],
                            destination,
                        )

                        download_rows.append({
                            "source_id": source["source_id"],
                            "kind": "ACS_CANDIDATE",
                            "original_url": link["url"],
                            "final_url": result["final_url"],
                            "status": result["download_status"],
                            "http_status": result["status_code"],
                            "content_type": result["content_type"],
                            "path": result["path"],
                            "bytes": result["bytes"],
                            "sha256": result["sha256"],
                            "error": result["error"],
                        })

                        time.sleep(SLEEP_SECONDS)

            else:
                download_rows.append({
                    "source_id": source["source_id"],
                    "kind": "LANDING_HTML",
                    "original_url": source["landing_url"],
                    "final_url": response.url,
                    "status": "HTTP_ERROR",
                    "http_status": response.status_code,
                    "content_type": response.headers.get(
                        "Content-Type",
                        "",
                    ),
                    "path": "",
                    "bytes": np.nan,
                    "sha256": "",
                    "error": f"HTTP {response.status_code}",
                })

        except requests.RequestException as error:
            download_rows.append({
                "source_id": source["source_id"],
                "kind": "LANDING_HTML",
                "original_url": source["landing_url"],
                "final_url": "",
                "status": "REQUEST_ERROR",
                "http_status": np.nan,
                "content_type": "",
                "path": "",
                "bytes": np.nan,
                "sha256": "",
                "error": str(error),
            })

        time.sleep(SLEEP_SECONDS)

    sources_csv = tables_dir / "rnr_primary_sources.csv"
    claims_csv = tables_dir / "rnr_source_claims.csv"
    links_csv = tables_dir / "rnr_candidate_links.csv"
    downloads_csv = tables_dir / "rnr_download_log.csv"

    pd.DataFrame(source_rows).to_csv(
        sources_csv,
        index=False,
    )

    pd.DataFrame(CLAIMS).to_csv(
        claims_csv,
        index=False,
    )

    pd.DataFrame(link_rows).to_csv(
        links_csv,
        index=False,
    )

    pd.DataFrame(download_rows).to_csv(
        downloads_csv,
        index=False,
    )

    raw_template = templates_dir / "neo_raw_samples_template.csv"
    surface_template = templates_dir / "neo_barrier_scenario_template.csv"
    request_file = templates_dir / "author_data_request.txt"

    make_raw_template(raw_template)
    make_surface_template(surface_template)
    make_author_request(request_file)

    readme = outdir / "README.md"

    readme.write_text(
        """# RNR NEO source archive

## Purpose

This directory is a provenance archive for the RNR QVSM/NEO model.

It does not contain invented QM/MM data.

## Files

- `tables/rnr_primary_sources.csv`: primary references and roles.
- `tables/rnr_source_claims.csv`: claims with scope and evidence state.
- `tables/rnr_candidate_links.csv`: candidate SI/data links found on landing pages.
- `tables/rnr_download_log.csv`: exact download attempts, status and SHA-256.
- `templates/neo_raw_samples_template.csv`: required schema for real configuration-level NEO reweighting data.
- `templates/neo_barrier_scenario_template.csv`: published scenario-level barrier/NQE shifts.
- `templates/author_data_request.txt`: request template for raw configuration-level data.

## Critical distinction

Per configuration/snapshot:

- `E_DFT_kcal_mol`
- `E_NEO_H_kcal_mol`
- `E_NEO_D_kcal_mol`
- `weight_classical`
- `R_OO_A`
- `delta_A`

Derived per configuration:

- `deltaE_NQE_H = E_NEO_H - E_DFT`
- `deltaE_NQE_D = E_NEO_D - E_DFT`

Derived from the ensemble/free-energy surface:

- `deltaG_classical_kcal_mol`
- `deltaG_quantum_H_kcal_mol`
- `deltaG_quantum_D_kcal_mol`
- `Q_deltaDeltaG_NQE_kcal_mol`
- `k_H`, `k_D`, `KIE`

Do not put a barrier into every snapshot row.
Do not assign `weight=1` unless a documented ensemble really has equal weights.
Do not infer D from an H-only NQE result.
""",
        encoding="utf-8",
    )

    files_for_manifest = [
        sources_csv,
        claims_csv,
        links_csv,
        downloads_csv,
        raw_template,
        surface_template,
        request_file,
        readme,
    ]

    file_hashes = {}

    for path in files_for_manifest:
        file_hashes[str(path.relative_to(outdir))] = sha256_file(path)

    manifest = {
        "project": "RNR-NEO-QM-MM-source-archive",
        "purpose": (
            "Source collection and data-template generation for "
            "reproducible QVSM-RNR NEO reweighting."
        ),
        "references": [
            source["doi"]
            for source in PRIMARY_SOURCES
        ],
        "NEO_reweighting_definition": {
            "deltaE_NQE_H": "E_NEO_H - E_DFT",
            "deltaE_NQE_D": "E_NEO_D - E_DFT",
            "w_H": "w_classical * exp[-beta*deltaE_NQE_H]",
            "w_D": "w_classical * exp[-beta*deltaE_NQE_D]",
            "G_H": "-RT ln(sum_bin(w_H)) + constant",
            "G_D": "-RT ln(sum_bin(w_D)) + constant",
        },
        "published_scenario_values": {
            "e4_alpha_E623_mediated_Q_H_kcal_mol": 0.8,
            "e4_alpha_prime_direct_Q_H_kcal_mol": -4.4,
            "configuration_counts": {
                "alpha": 50000,
                "alpha_prime": 85800,
            },
        },
        "epistemic_rules": [
            "Scenario-level published values are not raw configuration energies.",
            "Raw NEO energy tables must come from SI, repository, or authors.",
            "No KIE is calculated unless isotope-resolved H and D data are available.",
            "No QED water-coherence parameter is treated as measured.",
        ],
        "files_sha256": file_hashes,
    }

    manifest_path = outdir / "manifest.json"

    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("\n=== ARCHIVIO FONTI RNR/NEO CREATO ===")
    print(f"Directory: {outdir}")
    print("\nFile principali:")
    print(sources_csv)
    print(claims_csv)
    print(links_csv)
    print(downloads_csv)
    print(raw_template)
    print(surface_template)
    print(request_file)
    print(manifest_path)

    print("\nPer vedere candidate SI/data links:")
    print(
        "column -s, -t < "
        f"{links_csv} | less -S"
    )

    print("\nPer vedere esiti download:")
    print(
        "column -s, -t < "
        f"{downloads_csv} | less -S"
    )


if __name__ == "__main__":
    main()
