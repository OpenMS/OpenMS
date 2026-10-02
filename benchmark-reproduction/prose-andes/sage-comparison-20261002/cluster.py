import sys
import pyopenms as p
scan = int(sys.argv[1]); lo, hi = float(sys.argv[2]), float(sys.argv[3])
exp = p.MSExperiment(); p.MzMLFile().load("eclipse_tmtpro_10855_nocand.mzML", exp)
spec = next(s for s in exp if int(s.getNativeID().split("scan=")[-1]) == scan)
s = p.MSSpectrum(spec); p.Normalizer().filterSpectrum(s); s.sortByPosition()
print("raw peaks in window:")
for pk in s:
    if lo <= pk.getMZ() <= hi: print(f"   {pk.getMZ():10.4f} {pk.getIntensity():8.4f}")
a = p.MSSpectrum(s)
p.Deisotoper.deisotopeAndSingleCharge(a, 20.0, True, 1, 3, False, 3, 10, False, True, True)
arrs = {x.getName(): x for x in a.getIntegerDataArrays()}
print("after deisotoping (no conversion), arrays", list(arrs))
for i, pk in enumerate(a):
    if lo - 2 <= pk.getMZ() <= hi: print(f"   {pk.getMZ():10.4f} {pk.getIntensity():8.4f}  " + "  ".join(f"{k}={arrs[k][i]}" for k in arrs))
