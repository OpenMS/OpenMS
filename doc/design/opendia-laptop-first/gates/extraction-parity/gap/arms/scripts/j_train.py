"""Post-hoc: train D11FDJ = D11FDI with featJ in place of featI (frozen trainer, same seeds), then W and T."""
import os, subprocess
N = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113'
PY = '/home/singjc/.jcode/scratch/torchenv/.venv/bin/python'
src = open(f'{N}/scripts_gap/train_gap.py').read()
old = "FA, FG, FI = (np.load(f'{O}/{RAW}/feat{k}.npy').astype(np.float64) if any(k in a[3:] for a in ARMS) else None for k in 'AGI')"
assert old in src
src = src.replace(old, old + "\nFJ = np.load(f'{O}/{RAW}/featJ.npy').astype(np.float64) if any('J' in a[3:] for a in ARMS) else None")
old2 = "    for key, F_ in (('A', FA), ('G', FG), ('I', FI)):"
assert old2 in src
src = src.replace(old2, "    for key, F_ in (('A', FA), ('G', FG), ('I', FI), ('J', FJ)):")
p = f'{N}/scripts_gap/_train_J.py'; open(p, 'w').write(src)
for raw in ('raw01', 'raw02'):
    subprocess.run([PY, p, raw, 'D11FDJ'], check=True, env={**os.environ, 'TRAIN_THREADS': '10', 'OMP_NUM_THREADS': '10'})
    subprocess.run([PY, f'{N}/scripts/make_w.py', raw, 'D11FDJ', 'lin_D11FDJ', 'D11FDJW'], check=True)
    subprocess.run([PY, f'{N}/scripts/calib_tail.py', raw, 'D11FDJW', 'D11FDJWT'], check=True)
