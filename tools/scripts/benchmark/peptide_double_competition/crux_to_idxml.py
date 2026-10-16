#!/usr/bin/env python3
"""Convert one-hit Crux TSV results and Crux's paired peptide list to OpenMS inputs."""

import argparse
import csv
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


def read_pairs(path):
    pairs = {}
    decoys = set()
    missing = set()
    with path.open(encoding="utf-8", newline="") as source:
        rows = csv.DictReader(source, delimiter="\t")
        if not {"target", "decoy(s)"}.issubset(rows.fieldnames or []):
            raise ValueError("Expected Crux tide-index.peptides.txt target and decoy(s) columns")
        for row in rows:
            target = row["target"]
            decoy = row["decoy(s)"]
            if not decoy:
                missing.add(target)
                continue
            if target in pairs or decoy in decoys or target == decoy:
                raise ValueError(f"Crux peptide list is not one-to-one at {target} -> {decoy}")
            pairs[target] = decoy
            decoys.add(decoy)
    return pairs, missing


def convert(results_path, pairs_path, idxml_path, pair_path):
    pairs, missing = read_pairs(pairs_path)
    paired_decoys = set(pairs.values())
    counts = Counter()
    rows = []
    protein_labels = {}

    with results_path.open(encoding="utf-8", newline="") as source:
        results = csv.DictReader(source, delimiter="\t")
        required = {"scan", "charge", "retention time", "spectrum precursor m/z", "xcorr score", "sequence", "protein id", "target/decoy", "original target sequence"}
        if not required.issubset(results.fieldnames or []):
            raise ValueError("Unexpected Crux tide-search.txt columns")
        for row in results:
            sequence = row["sequence"]
            label = row["target/decoy"]
            original_target = row["original target sequence"]
            if label == "decoy" and pairs.get(original_target) != sequence:
                raise ValueError(f"Observed Crux decoy does not match the index pairing: {original_target} -> {sequence}")
            if label == "target" and original_target != sequence:
                raise ValueError(f"Observed Crux target has unexpected original sequence: {original_target} != {sequence}")
            if label == "target" and sequence not in pairs:
                counts["unpaired_target_psms"] += 1
                continue
            if label == "decoy" and sequence not in paired_decoys:
                counts["unpaired_decoy_psms"] += 1
                continue
            if label not in {"target", "decoy"}:
                counts["unknown_label_psms"] += 1
                continue
            protein_tokens = [item for item in row["protein id"].split(",") if item]
            if not protein_tokens:
                counts["missing_protein_psms"] += 1
                continue
            for accession in protein_tokens:
                old = protein_labels.setdefault(accession, label)
                if old != label:
                    raise ValueError(f"Protein has conflicting target/decoy labels: {accession}")
            rows.append((row, protein_tokens, label))
            counts[f"{label}_psms"] += 1

    root = ET.Element("IdXML", {"version": "1.5"})
    ET.SubElement(root, "SearchParameters", {
        "id": "SP_0", "db": "OpenMS UPS1/yeast/cRAP FASTA searched with Crux Tide",
        "db_version": "", "taxonomy": "", "mass_type": "monoisotopic", "charges": "",
        "enzyme": "trypsin", "missed_cleavages": "1", "precursor_peak_tolerance": "25",
        "precursor_peak_tolerance_ppm": "true", "peak_mass_tolerance": "0.02",
        "peak_mass_tolerance_ppm": "false",
    })
    run = ET.SubElement(root, "IdentificationRun", {
        "date": "2026-10-12T00:00:00", "search_engine": "Crux Tide",
        "search_engine_version": "4.3.2", "search_parameters_ref": "SP_0",
    })
    protein_id = ET.SubElement(run, "ProteinIdentification", {
        "score_type": "xcorr score", "higher_score_better": "true", "significance_threshold": "0",
    })
    protein_refs = {}
    for accession, label in protein_labels.items():
        ref = f"PH_{len(protein_refs)}"
        protein_refs[accession] = ref
        hit = ET.SubElement(protein_id, "ProteinHit", {
            "id": ref, "accession": accession, "score": "0", "sequence": "",
        })
        ET.SubElement(hit, "UserParam", {"type": "string", "name": "target_decoy", "value": label})

    for row, accessions, label in rows:
        peptide_id = ET.SubElement(run, "PeptideIdentification", {
            "score_type": "xcorr score", "higher_score_better": "true", "significance_threshold": "0",
            "RT": row["retention time"], "MZ": row["spectrum precursor m/z"],
            "spectrum_reference": f"scan={row['scan']}",
        })
        hit = ET.SubElement(peptide_id, "PeptideHit", {
            "score": row["xcorr score"], "sequence": row["sequence"], "charge": row["charge"],
            "protein_refs": " ".join(protein_refs[a] for a in accessions),
        })
        ET.SubElement(hit, "UserParam", {"type": "string", "name": "target_decoy", "value": label})

    ET.indent(root, space="  ")
    ET.ElementTree(root).write(idxml_path, encoding="UTF-8", xml_declaration=True)
    with pair_path.open("w", encoding="ascii", newline="\n") as output:
        output.write("target_sequence\tdecoy_sequence\n")
        for target, decoy in sorted(pairs.items()):
            output.write(f"{target}\t{decoy}\n")
    counts["paired_sequences"] = len(pairs)
    counts["target_sequences_without_decoys"] = len(missing)
    counts["input_psms"] = sum(counts[key] for key in ("target_psms", "decoy_psms", "unpaired_target_psms", "unpaired_decoy_psms", "unknown_label_psms", "missing_protein_psms"))
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("paired_peptides", type=Path)
    parser.add_argument("idxml", type=Path)
    parser.add_argument("pairs", type=Path)
    args = parser.parse_args()
    counts = convert(args.results, args.paired_peptides, args.idxml, args.pairs)
    for key, value in sorted(counts.items()):
        print(f"{key}={value}")
    print(f"idxml={args.idxml}")
    print(f"pairs={args.pairs}")


if __name__ == "__main__":
    main()
