"""PLAN_GAP recorded item: funnel recovery per arm (evaluation only). For each raw, the human DIA-NN q<=0.01 IDs whose
OpenDIA winner is DIA-NN's peak under the control, split into accepted / rejected at the N@FDP_c 1 % cut; per arm:
how many control-rejected become accepted, how many control-accepted are lost, and the arm's own funnel stage 4."""
import os, numpy as np, pandas as pd, subprocess, sys
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114'
SRC = open(B + '/funnel.py').read().split("# Why stage 1 loses")[0].split("print('\\nFunnel")[0]
for raw in ('raw01', 'raw02'):
    res = {}
    for arm in ('D11FDWT', 'D11FDAWT', 'D11FDGWT', 'D11FDIWT', 'D11FDAGIWT'):
        os.environ['RAW'] = raw
        g = {}
        exec(SRC.replace("np.load(f'{D}/scores_D11FDWT.npy')", f"np.load('{B}/run/{raw}/scores_{arm}.npy')"), g)
        acc = g['near_is_winner'] & g['accepted']
        res[arm] = (g['near_is_winner'], acc, g['n1'], len(g['accepted'] - set(g['ref'].precursor_id)))
    ctl_nw, ctl_acc = res['D11FDWT'][0], res['D11FDWT'][1]
    ctl_rej = ctl_nw - ctl_acc
    print(f'{raw}: control right-peak winners {len(ctl_nw)}, accepted {len(ctl_acc)}, rejected {len(ctl_rej)}')
    for arm, (nw, acc, n1, extra) in res.items():
        print(f'  {arm:11s} N@FDP_c1% {n1:5d} | DIA-NN-peak accepted {len(acc):5d} | recovered from control-rejected {len(acc & ctl_rej):4d} '
              f'| lost of control-accepted {len(ctl_acc - acc):4d} | right-peak winners {len(nw)} | accepted not DIA-NN {extra}')
