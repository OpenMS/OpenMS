"""ProSE on the #10364 suite over time, against Sage and ANDES (Percolator PSMs at q <= 0.01, seeds 1/42/137).

Per file: mean over seeds. Per group: mean over the group's files. "All 20 files": sum of the per-file means.
Columns:
  2026-09-30 reference: #10335's branch configuration frozen in #10364 (archived per-file means).
  develop f5ea2d04 (arm base), 3a47278 (+#10391, d_base), 9517361 (+#10394, d_dedup),
  2f30b40 (+#10397: d_fc on low-resolution files, hfx_A2 and astral_A2, d_fc_hr on the other high-resolution files),
  f6c680f (+#10398, #10399: re_auto; #10398 cannot change this suite, every spectrum has a charge).
  Sage 0.14.7 and ANDES auto: the frozen 2026-09-30 measurements (archived per-file means).
"""
import csv
import json
import re
import sys
from pathlib import Path

from evaluate import GROUPS, MANIFEST, group_of

ROOT = Path(__file__).resolve().parent
ARCHIVE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/user/OpenMS/benchmark-reproduction/prose-andes/ISSUE-10364-ARCHIVE-20261002.md")


def archived():
    """Per-file Sage, ProSE (2026-09-30 reference) and ANDES means from the archived per-file table."""
    out = {}
    for line in ARCHIVE.read_text().splitlines():
        m = re.match(r"\| `(\w+)` \| [\d /]+ \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \| \d+ \|$", line)
        if m:
            out[m[1]] = dict(sage=float(m[2]), ref=float(m[3]), andes=float(m[4]))
    return out


def archived_groups():
    """Group means of the frozen measurements as published (unrounded per-file values), by group name."""
    names = {"Velos CID": "Velos CID", "HF-X HCD": "HF-X HCD", "Astral HCD": "Astral HCD", "Lumos HCD LFQ": "Lumos HCD LFQ",
             "Lumos CID TMT": "Lumos CID TMT", "Exploris 480 TMTpro": "Exploris 480 TMTpro", "timsTOF HT": "timsTOF HT"}
    out = {}
    for line in ARCHIVE.read_text().splitlines():
        m = re.match(r"\| ([^|`]+?) \| \d \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \| [+-][\d.]+% \| [+-][\d.]+% \|$", line)
        if m and m[1] in names:
            out[names[m[1]]] = dict(ref=float(m[2]), andes=float(m[3]), sage=float(m[4]))
    return out


def mean(d, *arms):
    for a in arms:
        p = ROOT / "results" / d / a / "summary.json"
        if p.exists():
            return json.loads(p.read_text())["mean_psms"]
    raise FileNotFoundError(f"{d}: {arms}")


COLUMNS = [("ref", None), ("f5ea2d04", ("base",)), ("3a47278", ("d_base",)), ("9517361", ("d_dedup",)),
           ("2f30b40", ("d_fc", "d_fc_hr")), ("f6c680f", ("re_auto",)), ("sage", None), ("andes", None)]


def main():
    arch = archived()
    groups = archived_groups()
    assert len(arch) == 20 and len(groups) == 7, (len(arch), len(groups))
    files = [r["id"] for r in MANIFEST]
    per_file = {}
    for d in files:
        per_file[d] = {c: (arch[d][c] if arms is None else mean(d, *arms)) for c, arms in COLUMNS}
    with (ROOT / "history_per_file.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["dataset", "group"] + [c for c, _ in COLUMNS])
        for d in files:
            w.writerow([d, group_of(d)] + [round(per_file[d][c], 2) for c, _ in COLUMNS])
    rows = []
    for _, g in GROUPS:
        fs = [d for d in files if group_of(d) == g]
        v = {c: sum(per_file[d][c] for d in fs) / len(fs) for c, _ in COLUMNS}
        v.update(groups[g])  # the published group means; the archived per-file values are rounded
        rows.append((g, len(fs), v))
    rows.append(("All 20 files (sum)", 20, {c: sum(per_file[d][c] for d in files) for c, _ in COLUMNS}))
    lines = ["group\tn\t" + "\t".join(c for c, _ in COLUMNS) + "\tnow_vs_sage_pct\tnow_vs_andes_pct"]
    for g, n, v in rows:
        lines.append(f"{g}\t{n}\t" + "\t".join(f"{v[c]:.1f}" for c, _ in COLUMNS)
                     + f"\t{100 * (v['f6c680f'] / v['sage'] - 1):+.1f}\t{100 * (v['f6c680f'] / v['andes'] - 1):+.1f}")
    text = "\n".join(lines) + "\n"
    (ROOT / "history_per_group.tsv").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
