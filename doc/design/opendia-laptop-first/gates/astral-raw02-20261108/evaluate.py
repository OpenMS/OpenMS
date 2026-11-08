"""Entrapment evaluation for the single raw 02 run (astral-entrapment-raw02-20261108/GATE.md).

Labels come from the library only. A target precursor is entrapment if every protein mapped to its peptide has an
accession containing "_p_target". Usage: evaluate.py <crossfit.tsv> <library.pqp> [--library-only]
"""
import sqlite3
import sys

import pandas as pd

# 300223 target precursors in the PQP. The raw 01 gate wrote 151239 human; the library gives 151240 once the one
# shared precursor is counted as human (r 0.98508 vs 0.98509). Recorded in GATE.md before the run.
R_EXPECTED = (148983, 151240)


def library_labels(pqp):
    con = sqlite3.connect(pqp)
    rows = pd.read_sql(
        "select pr.ID as precursor_id, pr.DECOY as decoy, p.PROTEIN_ACCESSION as accession "
        "from PRECURSOR pr join PRECURSOR_PEPTIDE_MAPPING pp on pp.PRECURSOR_ID = pr.ID "
        "join PEPTIDE_PROTEIN_MAPPING pm on pm.PEPTIDE_ID = pp.PEPTIDE_ID join PROTEIN p on p.ID = pm.PROTEIN_ID",
        con)
    con.close()
    targets = rows[rows.decoy == 0].copy()
    targets["entrapment_protein"] = targets.accession.str.contains("_p_target", regex=False)
    per = targets.groupby("precursor_id").entrapment_protein.all()
    return per  # True = entrapment, False = human (shared counts as human)


def main():
    labels = library_labels(sys.argv[2])
    n_e, n_h = int(labels.sum()), int((~labels).sum())
    print(f"library target precursors: entrapment={n_e} human={n_h} (expected {R_EXPECTED[0]} / {R_EXPECTED[1]})")
    if (n_e, n_h) != R_EXPECTED:
        raise SystemExit("library labels do not reproduce the gate counts")
    r = n_e / n_h
    print(f"r = {r:.4f}")
    if "--library-only" in sys.argv:
        return
    c = pd.read_csv(sys.argv[1], sep="\t")
    winners = c[(c.precursor_best == 1) & (c.is_decoy == 0)].copy()
    winners["q_value"] = pd.to_numeric(winners.q_value)
    winners["entrapment"] = winners.precursor_id.map(labels)
    if winners.entrapment.isna().any():
        raise SystemExit("winner precursor missing from library labels")
    for q in (0.001, 0.01, 0.05):
        acc = winners[winners.q_value <= q]
        ne = int(acc.entrapment.sum())
        nh = int((~acc.entrapment.astype(bool)).sum())
        n = ne + nh
        fdp_c = ne * (1 + 1 / r) / n if n else float("nan")
        fdp_lb = ne / n if n else float("nan")
        tag = "  <- gate quantity" if q == 0.01 else ""
        print(f"q<={q}: accepted={n} N_h={nh} N_e={ne} FDP_combined={fdp_c:.4f} FDP_lower_bound={fdp_lb:.4f}{tag}")
    acc = winners[winners.q_value <= 0.01]
    n = len(acc)
    ne = int(acc.entrapment.sum())
    fdp_c = ne * (1 + 1 / r) / n if n else float("nan")
    if n == 0:
        print("GATE: FAIL (no accepted targets, FDP not measurable)")
    else:
        print(f"GATE: {'PASS' if fdp_c <= 0.015 else 'FAIL'} (FDP_combined {fdp_c:.4f} vs 0.015)")
    anchors = pd.read_csv(sys.argv[1].replace("crossfit/crossfit.tsv", "cal/anchors.tsv"), sep="\t")
    a = anchors.precursor_id.map(labels)
    print(f"calibration anchors: {len(anchors)} entrapment share={a.mean():.4f} (purity lower bound {1 - 2 * a.mean():.4f})")


if __name__ == "__main__":
    main()
