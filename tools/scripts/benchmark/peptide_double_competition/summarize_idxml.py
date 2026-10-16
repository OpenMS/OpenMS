#!/usr/bin/env python3
"""Summarize peptide q-values and UPS-positive yields from benchmark idXML."""

import sys
import xml.etree.ElementTree as ET


def summarize(path):
    root = ET.parse(path).getroot()
    proteins = {}
    for hit in root.findall(".//ProteinHit"):
        proteins[hit.attrib["id"]] = hit.attrib.get("accession", "")
    rows = []
    for peptide_hit in root.findall(".//PeptideHit"):
        meta = {item.attrib.get("name"): item.attrib.get("value")
                for item in peptide_hit.findall("UserParam")}
        label = meta.get("target_decoy", "unknown")
        q_value = meta.get("peptide q-value")
        if q_value is None:
            continue
        refs = peptide_hit.attrib.get("protein_refs", "").split()
        accessions = [proteins.get(ref, "") for ref in refs]
        ups = any("HUMAN_UPS" in accession for accession in accessions)
        rows.append((peptide_hit.attrib.get("sequence", ""), label, float(q_value), ups))
    print(f"file={path} annotated_hits={len(rows)}")
    for threshold in (0.01, 0.05, 0.10):
        accepted = [row for row in rows if row[2] <= threshold]
        unique_targets = {row[0] for row in accepted if row[1] == "target"}
        unique_decoys = {row[0] for row in accepted if row[1] == "decoy"}
        ups_targets = {row[0] for row in accepted if row[1] == "target" and row[3]}
        print(f"q<={threshold:g}: accepted_hits={len(accepted)} target_peptides={len(unique_targets)} "
              f"decoy_peptides={len(unique_decoys)} target_UPS_peptides={len(ups_targets)}")


if __name__ == "__main__":
    for filename in sys.argv[1:]:
        summarize(filename)
