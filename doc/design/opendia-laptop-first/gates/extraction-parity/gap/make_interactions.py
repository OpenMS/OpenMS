"""Probe: signal-level interaction features for the scorer (label-free). For each window, the evidence inputs that
separate weak real peptides from decoys at matched score (PRECHECK_JOINT section 2) are multiplied by a centred signal
level, so the linear readout can weight library agreement more where signal is weak. Signal level = log transition
support count (heuristics col 0) standardised later by the scorer. Output inter_<raw>.npy [n, 5] float32."""
import sys, numpy as np
G = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114'
raw = sys.argv[1]
H = np.load(f'{G}/armE_{raw}/heuristics.npy'); F = np.load(f'{G}/armE_{raw}/feat12.npy'); I = np.load(f'{G}/cpp/featI_{raw}.npy')
lvl = np.log1p(H[:, 0])  # transition_support
lvl = lvl - np.median(lvl)
cols = [F[:, 2], I[:, 1], H[:, 7], -np.log1p(H[:, 3]), H[:, 6]]  # f12_libcos, i_libcos_resid, boundary_library_correlation, -log rt error, boundary_coelution
out = np.stack([c * lvl for c in cols], 1).astype(np.float32)
np.save(f'/home/singjc/.jcode/scratch/e2e/joint/inter_{raw}.npy', out)
print(raw, out.shape, np.isfinite(out).all())
