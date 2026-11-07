#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
"""Sample raw-trace training windows with external seed labels, not engineered inputs.

This pilot requires the identical PQP used by OpenDIA and the scored teacher OSW.
Decoy TRAML_ID must be DECOY_<target TRAML_ID>. Group by the parent target's
unmodified sequence, keeping charge states, modifications and decoys together.
Teacher identifications only label windows. Their scores are never model inputs.
"""
import argparse
from contextlib import ExitStack
import csv
import hashlib
import math
from pathlib import Path
import sqlite3


def prepare(library, teacher, candidates, output, max_samples=6000, seed_q=0.001, all_windows=None):
    if max_samples < 6 or not 0 < seed_q <= 0.01:
        raise ValueError("Need at least six samples and a seed q threshold in (0, 0.01].")
    inputs = {Path(path).resolve() for path in (library, teacher, candidates)}
    if Path(output).resolve() in inputs or (all_windows and
            (Path(all_windows).resolve() in inputs or Path(all_windows).resolve() == Path(output).resolve())):
        raise ValueError("Window outputs must be distinct from inputs and each other.")
    with sqlite3.connect(f"file:{library}?mode=ro", uri=True) as connection:
        precursors = connection.execute(
            "SELECT p.ID,p.TRAML_ID,p.DECOY,pep.UNMODIFIED_SEQUENCE "
            "FROM PRECURSOR p JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID "
            "JOIN PEPTIDE pep ON pep.ID=m.PEPTIDE_ID"
        ).fetchall()
    targets = {name: sequence for _, name, decoy, sequence in precursors if not decoy}
    identity = {}
    for precursor_id, name, decoy, sequence in precursors:
        if decoy:
            if not name.startswith("DECOY_") or name[6:] not in targets:
                raise ValueError("Cannot resolve decoy parent. Supply explicit family grouping instead.")
            sequence = targets[name[6:]]
        value = (sequence, int(decoy), name)
        key = str(precursor_id)
        if not sequence or (key in identity and identity[key] != value):
            raise ValueError("Ambiguous or missing parent peptide identity.")
        identity[key] = value
    with sqlite3.connect(f"file:{teacher}?mode=ro", uri=True) as connection:
        # Reject identity remapping rather than silently assigning the wrong labels.
        for precursor_id, name, decoy in connection.execute("SELECT ID,TRAML_ID,DECOY FROM PRECURSOR"):
            if str(precursor_id) not in identity or identity[str(precursor_id)][1:] != (int(decoy), name):
                raise ValueError("Teacher OSW and PQP precursor IDs differ.")
        seeds = {}
        for precursor_id, left, right, apex in connection.execute(
            "SELECT f.PRECURSOR_ID,f.LEFT_WIDTH,f.RIGHT_WIDTH,f.EXP_RT FROM FEATURE f "
            "JOIN SCORE_MS2 s ON f.ID=s.FEATURE_ID JOIN PRECURSOR p ON p.ID=f.PRECURSOR_ID "
            "WHERE p.DECOY=0 AND s.RANK=1 AND s.QVALUE<=?", (seed_q,)
        ):
            if not all(math.isfinite(x) for x in (left, right, apex)) or left > right:
                raise ValueError("Non-finite or reversed teacher RT boundaries.")
            seeds.setdefault(str(precursor_id), []).append((left, right, apex))
    best = {}
    with ExitStack() as stack:
        handle = stack.enter_context(open(candidates, newline=""))
        all_writer = None
        if all_windows:
            all_writer = csv.writer(stack.enter_context(open(all_windows, "w", newline="")),
                                    delimiter="\t", lineterminator="\n")
            all_writer.writerow(["sample_id", "precursor_id", "group_id", "label", "rt_apex"])
        for row in csv.DictReader(handle, delimiter="\t"):
            precursor_id = row["precursor_id"]
            group, decoy, _ = identity[precursor_id]
            apex = float(row["rt_apex"])
            area = float(row["area"])
            if int(row["is_decoy"]) != decoy or not math.isfinite(apex) or not math.isfinite(area):
                raise ValueError("Invalid candidate coordinates or decoy status.")
            label = 0 if decoy else -1
            if decoy:
                rank = -area
            else:
                distances = [abs(apex - seed_apex) for left, right, seed_apex in seeds.get(precursor_id, [])
                             if left <= apex <= right]
                if distances:
                    label, rank = 1, min(distances)
            if all_writer:
                all_writer.writerow([row["candidate_row"], precursor_id, group, label, apex])
            if label == -1:
                continue  # An unconfirmed target is unknown, never a positive.
            value = (rank, row["candidate_row"], precursor_id, group, label, apex)
            if precursor_id not in best or value[:2] < best[precursor_id][:2]:
                best[precursor_id] = value
    pools = [[value for value in best.values() if value[4] == label] for label in (0, 1)]
    count = min(max_samples // 2, *(len(pool) for pool in pools))
    if count < 3:
        raise ValueError("Too few confirmed targets or decoys for training and held-out evaluation.")
    selected = []
    for pool in pools:
        pool.sort(key=lambda value: hashlib.sha256((value[2] + ":" + value[1]).encode()).digest())
        selected.extend(pool[:count])
    selected.sort(key=lambda value: (value[5], value[2]))
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "precursor_id", "group_id", "label", "rt_apex"])
        writer.writerows((row_id, precursor_id, group, label, apex)
                         for _, row_id, precursor_id, group, label, apex in selected)
    return {"positive_seeds": count, "decoys": count,
            "peptide_families": len({row[3] for row in selected}), "seed_q": seed_q}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", required=True)
    parser.add_argument("--teacher", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-samples", type=int, default=6000)
    parser.add_argument("--seed-q", type=float, default=0.001)
    parser.add_argument("--all-windows", help="Optional unsampled scoring window TSV, preserving unknown targets")
    print(prepare(**vars(parser.parse_args())))
