"""Replay ProSE's MS2 preprocessing (develop) on one spectrum and list which b/y ions of a peptide survive each step."""
import json, sys
import pyopenms as p
scan, pep = int(sys.argv[1]), sys.argv[2]
exp = p.MSExperiment(); p.MzMLFile().load("eclipse_tmtpro_10855_nocand.mzML", exp)
spec = next(s for s in exp if int(s.getNativeID().split("scan=")[-1]) == scan)
z = spec.getPrecursors()[0].getCharge()
seq = p.AASequence.fromString(pep)
ions = []
for kind, n in [("b", seq.size()), ("y", seq.size())]:
    for i in range(1, n):
        frag = seq.getPrefix(i) if kind == "b" else seq.getSuffix(i)
        for c in (1, 2):
            if c < z: ions.append((f"{kind}{i}{'+' * c}", frag.getMonoWeight(p.Residue.ResidueType.BIon if kind == "b" else p.Residue.ResidueType.YIon, c) / c))
def matched(s, tol=20):
    mz = [pk.getMZ() for pk in s]
    out = []
    for name, m in ions:
        if any(abs(x - m) / m * 1e6 <= tol for x in mz): out.append(name)
    return out
def show(label, s):
    m = matched(s)
    print(f"{label:28s} peaks {s.size():4d}  matched {len(m):2d}: {' '.join(m)}")
s = p.MSSpectrum(spec); show("raw", s)
tm = p.ThresholdMower(); prm = tm.getParameters(); prm.setValue("threshold", 0.05); tm.setParameters(prm); tm.filterSpectrum(s); show("ThresholdMower 0.05 (develop)", s)
p.Normalizer().filterSpectrum(s); s.sortByPosition()
d = p.MSSpectrum(s)
p.Deisotoper.deisotopeAndSingleCharge(d, 20.0, True, 1, 3, False, 3, 10, True)  # ProSEAlgorithm::preprocessSpectra_ call, library defaults after
show("deisotoped (1+ converted)", d)
a = p.MSSpectrum(s)
p.Deisotoper.deisotopeAndSingleCharge(a, 20.0, True, 1, 3, False, 3, 10, False, True)
ch = a.getIntegerDataArrays()[0] if a.getIntegerDataArrays() else None
names = dict(ions)
mz = [pk.getMZ() for pk in s]
for name, m in ions:
    hit = [x for x in mz if abs(x - m) / m * 1e6 <= 20]
    if not hit: continue
    j = min(range(a.size()), key=lambda k: abs(a[k].getMZ() - hit[0]))
    kept = abs(a[j].getMZ() - hit[0]) / hit[0] * 1e6 <= 20
    print(f"  {name:6s} {hit[0]:10.4f}  after deisotoping (no charge conversion): {'kept, charge ' + str(ch[j]) if kept and ch is not None else 'kept' if kept else 'REMOVED'}")
