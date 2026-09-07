#!/usr/bin/env python3
"""Acceptance-criteria analysis for bigbio/andes#64 on one fraction.

    offset_analysis.py <arm_dir> [<baseline_arm_dir>] [--testset issue64_testset.tsv] [--truth pglyco2.tsv.gz]

Per arm (from <arm>.glyco.pin + <arm>.t.psms):
  * accepted PSMs @1% by EFFECTIVE precursor offset = MonoShift (0 when the PIN has no
    Mono columns) + isotope_error — the number of isotopes between the recorded
    precursor and the searched monoisotope;
  * per offset tier: reference coverage (pGlyco2 / MSFragger scans), entrapment hits,
    HexNAc3 share and NeuGc>=1 share of the accepted compositions (the "chemically
    wrong tier" test), and the 84-scan test-set outcomes;
  * with a baseline arm: reference spectra gained / lost, and FLIPS — spectra
    CONFIRMED in the baseline that are no longer CONFIRMED here.
"""
import csv, gzip, re, sys, os
from collections import Counter, defaultdict

GLY = re.compile(r"HexNAc(\d+)Hex(\d+)Fuc(\d+)NeuAc(\d+)NeuGc(\d+)")

def read_pin(path):
    rows = {}
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        ix = {c: i for i, c in enumerate(hdr)}
        for r in rd:
            if len(r) < len(hdr) - 1:
                continue
            rows[r[ix["SpecId"]]] = r
    return hdr, ix, rows

def read_psms(path, q=0.01):
    acc = {}
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        qi = hdr.index("q-value"); pi = hdr.index("peptide")
        for r in rd:
            try:
                qv = float(r[qi])
            except (ValueError, IndexError):
                continue
            if qv <= q:
                acc[r[0]] = r[pi]
    return acc

def read_buckets(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["scan"])] = r["bucket"]
    return out

def scan_of(specid):
    m = re.search(r"scan=(\d+)", specid)
    return int(m.group(1)) if m else None

def analyse(arm_dir, testset, truth_scans, msf_scans):
    name = os.path.basename(arm_dir.rstrip("/"))
    hdr, ix, pin = read_pin(os.path.join(arm_dir, f"{name}.glyco.pin"))
    acc = read_psms(os.path.join(arm_dir, f"{name}.t.psms"))
    has_mono = "MonoShift" in ix
    by_off = Counter(); ref = Counter(); msf = Counter(); entrap = Counter()
    hexnac3 = Counter(); neugc = Counter(); shifted_rows = Counter()
    n_accepted = 0
    for sid in acc:
        r = pin.get(sid)
        if r is None:
            continue
        n_accepted += 1
        iso = int(float(r[ix["isotope_error"]]))
        shift = int(float(r[ix["MonoShift"]])) if has_mono else 0
        k = shift + iso
        by_off[k] += 1
        if shift:
            shifted_rows[shift] += 1
        scan = scan_of(sid)
        if scan in truth_scans: ref[k] += 1
        if scan in msf_scans: msf[k] += 1
        if "ENTRAP_" in r[-1]: entrap[k] += 1
        m = GLY.search(r[ix["Peptide"]])
        if m:
            n, h, f, a, g = (int(x) for x in m.groups())
            if n == 3: hexnac3[k] += 1
            if g >= 1: neugc[k] += 1
    print(f"\n=== {name}: {n_accepted} accepted glycoPSMs @1% (Mono columns: {has_mono})")
    if has_mono:
        print(f"    rows with an applied MonoShift among accepted: {dict(sorted(shifted_rows.items()))}")
    print("offset\taccepted\tin_pGlyco2\tin_MSFragger\tentrap\tHexNAc3%\tNeuGc>=1%")
    for k in sorted(by_off):
        n = by_off[k]
        print(f"{k:+d}\t{n}\t{ref[k]}\t{msf[k]}\t{entrap[k]}\t{100*hexnac3[k]/n:.1f}\t{100*neugc[k]/n:.1f}")
    buckets = read_buckets(os.path.join(arm_dir, "buckets_pglyco2.tsv"))
    if testset:
        c = Counter(buckets.get(s, "?") for s in testset)
        by_issue_offset = defaultdict(Counter)
        for s, row in testset.items():
            by_issue_offset[round(float(row["offset_Da"]))][buckets.get(s, "?")] += 1
        print(f"84-scan test set: {dict(c)}")
        for k in sorted(by_issue_offset):
            print(f"   issue offset +{k}: {dict(by_issue_offset[k])}")
    return buckets, acc, pin, ix

def main():
    args = sys.argv[1:]
    testset_p = None; truth_p = None
    if "--testset" in args:
        i = args.index("--testset"); testset_p = args[i+1]; del args[i:i+2]
    if "--truth" in args:
        i = args.index("--truth"); truth_p = args[i+1]; del args[i:i+2]
    arm = args[0]; base = args[1] if len(args) > 1 else None
    testset = {}
    if testset_p:
        with open(testset_p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                testset[int(r["scan"])] = r
    truth_scans = set(); msf_scans = set()
    G = "/home/user/andes/docs/benchmarks/glyco/truth/"
    for path, dest in ((G + "pglyco2_mouse_liver.tsv.gz", truth_scans), (G + "msfragger_mouse_liver.tsv.gz", msf_scans)):
        with gzip.open(path, "rt") as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                if r["run"] == "MouseLiver-Z-T-1":
                    dest.add(int(r["scan"]))
    b_arm, acc_arm, pin_arm, ix_arm = analyse(arm, testset, truth_scans, msf_scans)
    if base:
        b_base, acc_base, pin_base, ix_base = analyse(base, testset, truth_scans, msf_scans)
        conf_a = {s for s, b in b_arm.items() if b == "CONFIRMED"}
        conf_b = {s for s, b in b_base.items() if b == "CONFIRMED"}
        print(f"\n=== {os.path.basename(arm)} vs {os.path.basename(base)} (pGlyco2 reference)")
        print(f"confirmed: {len(conf_b)} -> {len(conf_a)}; gained {len(conf_a - conf_b)}, lost (flips) {len(conf_b - conf_a)}")
        flips = sorted(conf_b - conf_a)
        if flips:
            print("flipped scans:", ", ".join(f"{s}:{b_arm.get(s)}" for s in flips[:40]), "..." if len(flips) > 40 else "")
        # accepted-set overlap
        sa = set(scan_of(s) for s in acc_arm); sb = set(scan_of(s) for s in acc_base)
        print(f"accepted scans: base {len(sb)}, arm {len(sa)}, both {len(sa & sb)}, only arm {len(sa - sb)}, only base {len(sb - sa)}")
        # same-scan peptidoform agreement between the arms on shared accepted scans
        same = sum(1 for s in acc_arm if s in acc_base and acc_arm[s] == acc_base[s])
        shared = sum(1 for s in acc_arm if s in acc_base)
        print(f"same-scan peptidoform agreement between arms: {same}/{shared} = {100*same/max(shared,1):.1f}%")

if __name__ == "__main__":
    main()
