"""A+B fix, step 1 (label-free): extended RT calibration. Usage: extend_cal.py <in cal.json> <out cal.json>.
OpenDIA maps library RT through the calibration knots piecewise-linearly and holds the end knot value constant outside
them (RTMap::operator()). The knots cover library RT 0.098..0.752 (raw 01), but the library spans -0.049..0.940, so 44 % of
library precursors get the end value. This adds one knot below and one above, at the library RT range ends
(-0.049 and 0.940, the PQP min / max), placed on the line through the end knot and the 3rd knot from that end. Inside
the original knots the map is unchanged, and the RT tolerance, ms2 and im blocks are unchanged. Decoys share their
target's library RT, so the change is the same for targets and decoys."""
import json, sys
src, dst = sys.argv[1:3]
LO, HI, K = -0.048912644386292, 0.940189480781555, 3
c = json.load(open(src)); rt = c['rt']; x, y = list(rt['knots_library_rt']), list(rt['knots_run_rt_s'])
sl = (y[K] - y[0]) / (x[K] - x[0]); sh = (y[-1] - y[-1 - K]) / (x[-1] - x[-1 - K])
x = [LO] + x + [HI]; y = [y[0] + sl * (LO - x[1])] + y + [y[-1] + sh * (HI - x[-2])]
assert all(b > a for a, b in zip(x, x[1:])) and all(b >= a for a, b in zip(y, y[1:]))
rt['knots_library_rt'], rt['knots_run_rt_s'] = x, y
c['rt_extension'] = {'rule': 'linear end extrapolation through end knot and 3rd knot from the end', 'library_rt_range': [LO, HI],
                     'low_slope_s': sl, 'high_slope_s': sh}
json.dump(c, open(dst, 'w'), indent=1)
print(f'{dst}: knots {len(x)}, added ({LO:.3f}, {y[0]:.1f} s) and ({HI:.3f}, {y[-1]:.1f} s); slopes {sl:.0f} / {sh:.0f} s per unit')
