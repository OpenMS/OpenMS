"""C++ port check (post-hoc, descriptive): D11FDI retrained with featI from the C++ OpenDIAInterference tool (featC.npy,
same 4 columns, same row order) in place of the Python featI, at the frozen seed and seed offsets 1-4, against the
own-seed D11FD control already computed (scores_D11FDWT / scores_D11FD_s<o>WT). Usage: cpp_check.py <offset>."""
import os, subprocess, sys
off = int(sys.argv[1])
N = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113'
PY = '/home/singjc/.jcode/scratch/torchenv/.venv/bin/python'
src = open(f'{N}/scripts_gap/train_gap.py').read()
old = "FA, FG, FI = (np.load(f'{O}/{RAW}/feat{k}.npy').astype(np.float64) if any(k in a[3:] for a in ARMS) else None for k in 'AGI')"
assert old in src
src = src.replace(old, old + "\nFC = np.load(f'{O}/{RAW}/featC.npy').astype(np.float64) if any('C' in a[4:] for a in ARMS) else None")
old2 = "    for key, F_ in (('A', FA), ('G', FG), ('I', FI)):"
assert old2 in src
src = src.replace(old2, "    for key, F_ in (('A', FA), ('G', FG), ('I', FI), ('C', FC)):")
# 'C' is also the featPC flag letter in train_v4; arm D11FDC: only our block (FPC needs 'C' in a[3:] -> guard below).
src = src.replace("FPC = np.load(f'{O}/{RAW}/featPC.npy').astype(np.float64) if any(('P' in a[3:]) or ('C' in a[3:]) for a in ARMS) else None", "FPC = None")
src = src.replace("    for key, cols in PC_COLS.items():\n        if key in arm[3:]: Hx = np.c_[Hx, FPC[:, cols]]\n", "")
if off:
    src = src.replace('100 + k, valid', f'100 + k + {off}, valid').replace('200 + 10 * k + r, valid', f'200 + 10 * k + r + {off}, valid')
tag = f'_s{off}' if off else ''
src = src.replace("np.save(f'{O}/{RAW}/scores_{arm}.npy', scores); np.save(f'{O}/{RAW}/lin_{arm}.npy', lin)",
                  f"np.save(f'{{O}}/{{RAW}}/scores_{{arm}}{tag}.npy', scores); np.save(f'{{O}}/{{RAW}}/lin_{{arm}}{tag}.npy', lin)")
assert 'FPC = None' in src
p = f'{N}/scripts_gap/_train_C{off}.py'; open(p, 'w').write(src)
for raw in ('raw01', 'raw02'):
    subprocess.run([PY, p, raw, 'D11FDC'], check=True, env={**os.environ, 'TRAIN_THREADS': '6', 'OMP_NUM_THREADS': '6'})
    subprocess.run([PY, f'{N}/scripts/make_w.py', raw, f'D11FDC{tag}', f'lin_D11FDC{tag}', f'D11FDC{tag}W'], check=True)
    subprocess.run([PY, f'{N}/scripts/calib_tail.py', raw, f'D11FDC{tag}W', f'D11FDC{tag}WT'], check=True)
