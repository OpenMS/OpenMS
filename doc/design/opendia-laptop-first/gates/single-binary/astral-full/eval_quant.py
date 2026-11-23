#!/usr/bin/env python3
"""Descriptive OpenDIA area repeatability. Standard library only; no imputation."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics as st
import sqlite3


def read_run(path, cutoff):
    values, seen, counts = {}, set(), dict(rows=0, decoys=0, rejected_q=0, nonpositive=0, nonfinite=0)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"run", "precursor_id", "is_decoy", "sequence", "charge", "area", "q_value"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f"{path}: expected single-run precursor TSV columns {sorted(required)}")
        runs = set()
        for row in reader:
            counts["rows"] += 1
            runs.add(row["run"])
            if row["is_decoy"] not in {"0", "1"}:
                raise ValueError("is_decoy must be 0 or 1")
            if row["is_decoy"] == "1":
                counts["decoys"] += 1
                continue
            key = row["precursor_id"]
            if not key or not row["sequence"] or int(row["charge"]) <= 0 or key in seen:
                raise ValueError(f"{path}: invalid or duplicate target precursor {key}")
            seen.add(key)
            q, area = float(row["q_value"]), float(row["area"])
            if not math.isfinite(q) or not 0 <= q <= 1:
                raise ValueError(f"{path}: invalid q_value")
            if q > cutoff:
                counts["rejected_q"] += 1
                continue
            if not math.isfinite(area):
                counts["nonfinite"] += 1
                area = None
            elif area <= 0:
                counts["nonpositive"] += 1
                area = None
            values[key] = (row["sequence"], int(row["charge"]), area)
        if len(runs) > 1:
            raise ValueError(f"{path}: multiple raw runs in single-run input")
    return values, counts


def correlation(xs, ys):
    if len(xs) < 2 or len(set(xs)) < 2 or len(set(ys)) < 2:
        return None
    return st.correlation(xs, ys)


def describe(matrix, labels):
    keys = sorted(set().union(*(set(x) for x in matrix)))
    common = [k for k in keys if all(x.get(k) is not None for x in matrix)]
    medians = [st.median([math.log2(x[k]) for k in common]) if common else None for x in matrix]
    center = st.median(medians) if common else None
    # Subtract median log-area on the all-run complete set, retain original scale center.
    offsets = [m - center if m is not None else None for m in medians]
    rows = []
    for key in keys:
        raw = [x.get(key) for x in matrix]
        normalized = [2 ** (math.log2(v) - o) if v is not None and o is not None else None
                      for v, o in zip(raw, offsets)]
        complete = all(v is not None for v in raw)
        cv = lambda vs: st.stdev(vs) / st.mean(vs) if complete and len(vs) > 1 else None
        observed = [v for v in raw if v is not None]
        rows.append(dict(feature=key, detected=sum(key in x for x in matrix), positive=len(observed),
                         raw=raw, median_normalized=normalized, cv_raw=cv(raw),
                         cv_median_normalized=cv(normalized) if common else None,
                         observed_median_log2=st.median(map(math.log2, observed)) if observed else None))
    pairs = []
    for i, j in itertools.combinations(range(len(matrix)), 2):
        a, b = set(matrix[i]), set(matrix[j])
        shared = sorted(k for k in a & b if matrix[i][k] is not None and matrix[j][k] is not None)
        pairs.append(dict(a=labels[i], b=labels[j], intersection=len(a & b), union=len(a | b),
                          jaccard=len(a & b) / len(a | b) if a | b else None,
                          positive_pairs=len(shared), log2_pearson=correlation(
                              [math.log2(matrix[i][k]) for k in shared],
                              [math.log2(matrix[j][k]) for k in shared])))
    return dict(labels=labels, union=len(keys), positive_complete=len(common),
                detection_completeness={str(n): sum(r["detected"] == n for r in rows)
                                        for n in range(1, len(labels) + 1)},
                positive_completeness={str(n): sum(r["positive"] == n for r in rows)
                                      for n in range(len(labels) + 1)},
                normalization_log2_offsets=offsets, pairs=pairs, features=rows)


def population_labels(path):
    """Post-search descriptive strata only, never score or q-value calibration."""
    with sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True) as db:
        rows = db.execute("""SELECT p.ID, max(instr(pr.PROTEIN_ACCESSION, '_p_target')>0)
            FROM PRECURSOR p LEFT JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID
            LEFT JOIN PEPTIDE_PROTEIN_MAPPING pp ON pp.PEPTIDE_ID=m.PEPTIDE_ID
            LEFT JOIN PROTEIN pr ON pr.ID=pp.PROTEIN_ID WHERE p.DECOY=0 GROUP BY p.ID""").fetchall()
    if any(entrapment is None for _, entrapment in rows):
        raise ValueError("population grading requires protein mapping for every library target")
    return {str(pid): bool(entrapment) for pid, entrapment in rows}


def evaluate(manifest, cutoff, library=None, population="all"):
    if population not in {"all", "human"} or population == "human" and library is None:
        raise ValueError("human population requires --library; population must be all or human")
    labels = population_labels(library) if library is not None else None
    with open(manifest, newline="") as handle:
        entries = list(csv.DictReader(handle, delimiter="\t"))
    if not entries or set(entries[0]) != {"raw", "seed", "path"}:
        raise ValueError("manifest requires exactly raw, seed, path columns")
    data, audit, metadata = {}, [], {}
    for entry in entries:
        raw, seed = entry["raw"], entry["seed"]
        if not raw or not seed or (raw, seed) in data:
            raise ValueError("empty or duplicate raw/seed")
        path = Path(manifest).resolve().parent / entry["path"]
        values, counts = read_run(path, cutoff)
        if labels is not None and any(key not in labels for key in values):
            raise ValueError("accepted precursor ID not a target in grading library")
        counts["population_excluded"] = sum(labels[key] for key in values) if population == "human" else 0
        if population == "human":
            values = {key: value for key, value in values.items() if not labels[key]}
        for key, (sequence, charge, _) in values.items():
            if key in metadata and metadata[key] != (sequence, charge):
                raise ValueError(f"inconsistent library identity for {key}")
            metadata[key] = (sequence, charge)
        precursor = {k: v[2] for k, v in values.items()}
        peptide = {}
        for sequence, _, area in values.values():
            peptide.setdefault(sequence, [])
            if area is not None:
                peptide[sequence].append(area)
        peptide = {k: math.fsum(v) if v else None for k, v in peptide.items()}
        data[raw, seed] = dict(precursor=precursor, peptide=peptide)
        audit.append(dict(raw=raw, seed=seed, path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(), **counts))
    raws = sorted({r for r, _ in data})
    seeds = sorted({s for _, s in data})
    if len(raws) != 3 or set(data) != set(itertools.product(raws, seeds)):
        raise ValueError("require exactly three raws and the same seed set for every raw")
    output = dict(q_cutoff=cutoff, population=population,
                  population_definition="all accepted targets" if population == "all" else
                  "accepted targets without any protein accession containing literal _p_target (human proxy, may include contaminants)",
                  grading_library=str(Path(library).resolve()) if library is not None else None,
                  audit=audit, across_raws={}, within_raw_seeds={})
    for level in ("precursor", "peptide"):
        output["across_raws"][level] = {s: describe([data[r, s][level] for r in raws], raws) for s in seeds}
        output["within_raw_seeds"][level] = {r: describe([data[r, s][level] for s in seeds], seeds) for r in raws}
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", help="TSV: raw, seed, path (relative to manifest)")
    parser.add_argument("--q-cutoff", type=float, default=0.01)
    parser.add_argument("--library", help="frozen search PQP for post-search population grading")
    parser.add_argument("--population", choices=("all", "human"), default="all")
    args = parser.parse_args()
    if not math.isfinite(args.q_cutoff) or not 0 <= args.q_cutoff <= 1:
        parser.error("q cutoff must be in [0, 1]")
    print(json.dumps(evaluate(args.manifest, args.q_cutoff, args.library, args.population), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
