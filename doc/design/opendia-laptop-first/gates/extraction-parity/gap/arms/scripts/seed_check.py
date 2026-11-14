"""Post-hoc seed sensitivity (descriptive, not part of PLAN_GAP's ruling). train_gap.py with every training seed shifted
by SEED_OFFSET (offset 0 = the frozen run), for D11FD and D11FDI, then make_w and calib_tail as the chain, written as
scores_<arm>_s<offset>WT. Usage (inside the SIF, run root bound as the chain): seed_check.py <offset>."""
import os, subprocess, sys
off = int(sys.argv[1])
N = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113'
PY = '/home/singjc/.jcode/scratch/torchenv/.venv/bin/python'
src = open(f'{N}/scripts_gap/train_gap.py').read()
src = src.replace('100 + k, valid', f'100 + k + {off}, valid').replace('200 + 10 * k + r, valid', f'200 + 10 * k + r + {off}, valid')
src = src.replace("np.save(f'{O}/{RAW}/scores_{arm}.npy', scores); np.save(f'{O}/{RAW}/lin_{arm}.npy', lin)",
                  f"np.save(f'{{O}}/{{RAW}}/scores_{{arm}}_s{off}.npy', scores); np.save(f'{{O}}/{{RAW}}/lin_{{arm}}_s{off}.npy', lin)")
assert f'_s{off}.npy' in src and f'+ {off}, valid' in src
path = f'{N}/scripts_gap/_train_seed{off}.py'; open(path, 'w').write(src)
for raw in ('raw01', 'raw02'):
    subprocess.run([PY, path, raw, 'D11FD', 'D11FDI'], check=True, env={**os.environ, 'TRAIN_THREADS': '6', 'OMP_NUM_THREADS': '6'})
    for a in ('D11FD', 'D11FDI'):
        subprocess.run([PY, f'{N}/scripts/make_w.py', raw, f'{a}_s{off}', f'lin_{a}_s{off}', f'{a}_s{off}W'], check=True)
        subprocess.run([PY, f'{N}/scripts/calib_tail.py', raw, f'{a}_s{off}W', f'{a}_s{off}WT'], check=True)
