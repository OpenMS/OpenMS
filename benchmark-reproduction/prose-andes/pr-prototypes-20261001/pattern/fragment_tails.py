"""Fragment mass-error tails per file: matched ions of accepted target PSMs vs decoy PSMs.

Reads the develop (base) search idXML (fragment_annotation: observed m/z, relative intensity, charge, ion name of
every matched peak in the preprocessed spectrum) and Percolator's seed-1 target (q <= 0.01) and decoy tables.
Theoretical m/z come from TheoreticalSpectrumGenerator. For the calibrated fragment window w of #10379's auto
mode (c_default), reports the share of matched ions beyond w for targets (true ions lost) and decoys (random
matches removed), overall and for weak (< 5% of base peak) and strong (>= 20%) peaks.
"""
import csv, gzip, json, re, statistics, sys
from pathlib import Path
import pyopenms as oms

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = json.loads((ROOT.parent / "prose-andes-reproduction/broad_benchmark/manifest.json").read_text())


def scan_of(ref):
    m = re.search(r"scan=(\d+)", ref)
    return int(m.group(1)) if m else None


def stripped(pep):
    plain = re.sub(r"\[[^\]]*\]", "", pep)
    core = plain.split(".")[1] if plain.count(".") == 2 else plain
    return re.sub(r"[^A-Z]", "", core)


def accepted(path, qmax):
    out = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if float(row["q-value"]) <= qmax:
                sid = row["PSMId"].rsplit("_", 2)
                out.append((int(sid[-2]), stripped(row["peptide"])))
    return out


def main(datasets):
    tsg = oms.TheoreticalSpectrumGenerator()
    p = tsg.getParameters(); p.setValue("add_metainfo", "true"); tsg.setParameters(p)
    rows = []
    for d in datasets:
        base = ROOT / "results" / d / "base"
        tmp = ROOT / "pattern" / f"{d}.idXML"
        if not tmp.exists():
            tmp.write_bytes(gzip.open(base / "search.idXML.gz").read())
        prots, peps = [], oms.PeptideIdentificationList() if hasattr(oms, "PeptideIdentificationList") else []
        oms.IdXMLFile().load(str(tmp), prots, peps)
        by_scan = {}
        for pid in peps:
            s = scan_of(pid.getMetaValue("spectrum_reference") if pid.metaValueExists("spectrum_reference") else pid.getSpectrumReference())
            by_scan[s] = pid.getHits()
        w = json.loads((ROOT / "results" / d / "c_default" / "summary.json").read_text())["resolved"]["fragment_mass_tolerance"]
        res = {}
        for label, path in [("target", base / "s1" / "target.tsv"), ("decoy", base / "s1" / "decoy.tsv")]:
            psms = accepted(path, 0.01 if label == "target" else 1.0)
            errs = []  # (ppm, rel_int)
            for scan, pep in psms:
                hit = next((h for h in by_scan.get(scan, []) if h.getSequence().toUnmodifiedString() == pep), None)
                if hit is None or not hit.getPeakAnnotations():
                    continue
                theo = oms.MSSpectrum()
                tsg.getSpectrum(theo, hit.getSequence(), 1, 3)
                names = [n.decode() if isinstance(n, bytes) else n for n in theo.getStringDataArrays()[0]]
                theo_mz = {n: theo[i].getMZ() for i, n in enumerate(names)}
                for ann in hit.getPeakAnnotations():
                    name = ann.annotation.decode() if isinstance(ann.annotation, bytes) else ann.annotation
                    if name not in theo_mz:
                        continue
                    errs.append(((ann.mz - theo_mz[name]) / theo_mz[name] * 1e6, ann.intensity))
            res[label] = errs
        t, dcy = res["target"], res["decoy"]
        if not t or not dcy:
            print(d, "no data", len(t), len(dcy), file=sys.stderr); continue
        te = [e for e, _ in t]
        med = statistics.median(te)
        sd = 1.4826 * statistics.median(abs(e - med) for e in te)
        frac = lambda xs: sum(abs(e) > w for e, _ in xs) / len(xs)
        weak = [x for x in t if x[1] < 0.05]; strong = [x for x in t if x[1] >= 0.2]
        rows.append(dict(dataset=d, window_ppm=round(w, 2), target_ions=len(t), decoy_ions=len(dcy),
                         median_ppm=round(med, 2), robust_sd_ppm=round(sd, 2),
                         sd_weak=round(1.4826 * statistics.median(abs(e - med) for e, _ in weak), 2) if weak else "",
                         sd_strong=round(1.4826 * statistics.median(abs(e - med) for e, _ in strong), 2) if strong else "",
                         target_beyond_w_pct=round(100 * frac(t), 2), weak_beyond_w_pct=round(100 * frac(weak), 2) if weak else "",
                         decoy_beyond_w_pct=round(100 * frac(dcy), 2),
                         weak_share_pct=round(100 * len(weak) / len(t), 1)))
        print(rows[-1], flush=True)
    with open(ROOT / "pattern" / "fragment_tails.tsv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main(sys.argv[1:])
