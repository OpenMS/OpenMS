"""Identity checks for the 2026-10-02 develop+fix arms.

- d_base (develop 3a47278 ProSE sources) vs ds1 (deisotope fix on f5ea2d04) on high-resolution files and
  vs base (f5ea2d04) on low-resolution files, which are not deisotoped: identical native TSV and PIN.
- d_79 (#10379 default) on the two low-resolution check files: identical to d_base (auto does not calibrate).
- d_priors (#10378): identical native TSV, and identical PIN after removing the ion_prior_* columns.
"""
import csv
import gzip
import hashlib
import io
from pathlib import Path

from evaluate import MANIFEST

ROOT = Path(__file__).resolve().parent


def text(p):
    return gzip.decompress(p.read_bytes()).decode() if p.suffix == ".gz" else p.read_text()


def pin_without(s, prefix):
    rows = list(csv.reader(io.StringIO(s), delimiter="\t"))
    keep = [i for i, n in enumerate(rows[0]) if not n.startswith(prefix)]
    out = io.StringIO()
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    for r in rows:
        w.writerow([r[i] for i in keep if i < len(r)])
    return out.getvalue()


def files(d, a):
    base = ROOT / "results" / d / a
    return {k: next((base / n for n in names if (base / n).exists()), None)
            for k, names in (("tsv", ("native.tsv.gz", "native.tsv")), ("pin", ("native.pin.gz", "input.pin")))}


def same(d, a, b, drop=None):
    fa, fb = files(d, a), files(d, b)
    if not all(fa.values()) or not all(fb.values()):
        return "missing"
    res = []
    for k in ("tsv", "pin"):
        ta, tb = text(fa[k]), text(fb[k])
        if k == "pin" and drop:
            ta, tb = pin_without(ta, drop), pin_without(tb, drop)
        res.append(f"{k} {'same' if hashlib.sha256(ta.encode()).digest() == hashlib.sha256(tb.encode()).digest() else 'DIFF'}")
    return ", ".join(res)


def main():
    lines = []
    for r in MANIFEST:
        d = r["id"]
        ref = "ds1" if r["fragment_unit"] == "ppm" else "base"
        lines.append(f"d_base vs {ref:4s} {d:22s} {same(d, 'd_base', ref)}")
    for d in ("velos_125_R1", "lumos_tmt_5058"):
        lines.append(f"d_79 vs d_base  {d:22s} {same(d, 'd_79', 'd_base')}")
    for r in MANIFEST:
        d = r["id"]
        lines.append(f"d_priors vs d_base (no ion_prior_*) {d:22s} {same(d, 'd_priors', 'd_base', drop='ion_prior_')}")
    (ROOT / "checks_d.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
