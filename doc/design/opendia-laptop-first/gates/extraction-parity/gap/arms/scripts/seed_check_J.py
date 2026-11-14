"""Post-hoc: seed sensitivity of variant J (block I without same-composition claimants), as seed_check.py for D11FDI.
Uses _train_J.py (written by j_train.py) with every training seed shifted by the offset; writes scores_D11FDJ_s<o>WT.
Usage (inside the SIF, run root bound as the chain): seed_check_J.py <offset>."""
import os, subprocess, sys
off = int(sys.argv[1])
N = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113'
PY = '/home/singjc/.jcode/scratch/torchenv/.venv/bin/python'
src = open(f'{N}/scripts_gap/_train_J.py').read()
src = src.replace('100 + k, valid', f'100 + k + {off}, valid').replace('200 + 10 * k + r, valid', f'200 + 10 * k + r + {off}, valid')
a = "np.save(f'{O}/{RAW}/scores_{arm}.npy', scores); np.save(f'{O}/{RAW}/lin_{arm}.npy', lin)"
assert a in src and f'+ {off}, valid' in src
src = src.replace(a, f"np.save(f'{{O}}/{{RAW}}/scores_{{arm}}_s{off}.npy', scores); np.save(f'{{O}}/{{RAW}}/lin_{{arm}}_s{off}.npy', lin)")
p = f'{N}/scripts_gap/_train_J_seed{off}.py'; open(p, 'w').write(src)
for raw in ('raw01', 'raw02'):
    subprocess.run([PY, p, raw, 'D11FDJ'], check=True, env={**os.environ, 'TRAIN_THREADS': '6', 'OMP_NUM_THREADS': '6'})
    subprocess.run([PY, f'{N}/scripts/make_w.py', raw, f'D11FDJ_s{off}', f'lin_D11FDJ_s{off}', f'D11FDJ_s{off}W'], check=True)
    subprocess.run([PY, f'{N}/scripts/calib_tail.py', raw, f'D11FDJ_s{off}W', f'D11FDJ_s{off}WT'], check=True)
