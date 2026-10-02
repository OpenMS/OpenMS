"""Emulate the TOPP_ProSE_DDA* tests (Bruker HeLa DDA, native 1% PSM FDR, PeptideHit floors) for two builds.

ProSE runs through the harness driver with the tests' Search:* overrides on top of the algorithm
defaults; FalseDiscoveryRate (PSM 1%, no protein FDR, decoys removed) then runs through pyOpenMS,
and the PeptideHits are counted as check_psm_floor.cmake does.
"""
import json, os, subprocess, sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
HERE = Path(__file__).resolve().parent
PU = HERE.parent / "pr-untangled"
BUNDLE = HERE.parent / "prose-andes-reproduction"
ENV = dict(os.environ, OMP_NUM_THREADS="2", OPENMS_DATA_PATH=str(BUNDLE / "nightly/pyopenms/share/OpenMS"),
           PYTHONPATH=str(BUNDLE / "nightly") + os.pathsep + str(BUNDLE / "deps"))
D = next(HERE.glob("*.d"))
BASE = {"precursor:mass_tolerance_lower": 20.0, "precursor:mass_tolerance_upper": 20.0, "precursor:mass_tolerance_unit": "ppm",
        "fragment:mass_tolerance": 20.0, "fragment:mass_tolerance_unit": "ppm", "enzyme": "Trypsin/P",
        "peptide:missed_cleavages": 2, "modifications:variable": ["Oxidation (M)", "Acetyl (Protein N-term)"], "decoys": "auto"}
CONFIGS = {"DDA": {}, "DDA_calibrated": {"calibration:enabled": "true"},
           "DDA_chunked_calibrated": {"calibration:enabled": "true", "database:chunk_size": 5000}}
FLOORS = {"DDA": 4350, "DDA_calibrated": 4650, "DDA_chunked_calibrated": None}

def run(build, cfg):
    out = HERE / "runs" / build / cfg
    out.mkdir(parents=True, exist_ok=True)
    (out / "params.json").write_text(json.dumps(dict(BASE, **CONFIGS[cfg]), indent=1))
    if not (out / "search.idXML").exists():
        with (out / "search.log").open("w") as log:
            subprocess.run([str(PU / "builds" / build / "prose"), str(D), str(HERE / "human_sp.fasta"), str(out / "params.json"),
                            str(out / "native.tsv"), str(out / "search.idXML")], env=ENV, stdout=log, stderr=log, check=True)
    code = f"""
import pyopenms as p
prots, peps = [], p.PeptideIdentificationList()
p.IdXMLFile().load({str(out / 'search.idXML')!r}, prots, peps)
fdr = p.FalseDiscoveryRate(); fp = fdr.getParameters(); fp.setValue("add_decoy_peptides", "false"); fdr.setParameters(fp)
fdr.apply(peps)
p.IDFilter.filterHitsByScore(peps, 0.01)
p.IDFilter.removeDecoyHits(peps)
print(sum(len(x.getHits()) for x in peps))
"""
    n = int(subprocess.run([sys.executable, "-c", code], env=ENV, capture_output=True, text=True, check=True).stdout.strip().splitlines()[-1])
    return build, cfg, n

if __name__ == "__main__":
    builds = sys.argv[1:]
    with ThreadPoolExecutor(2) as ex:
        res = list(ex.map(lambda a: run(*a), [(b, c) for c in CONFIGS for b in builds]))
    lines = []
    for c in CONFIGS:
        row = {b: n for b, cc, n in res if cc == c}
        base = row[builds[0]]
        lines.append(f"{c:24s} floor {FLOORS[c]}  " + "  ".join(f"{b} {row[b]} ({100 * (row[b] / base - 1):+.2f}%)" for b in builds))
    (HERE / "hela_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
