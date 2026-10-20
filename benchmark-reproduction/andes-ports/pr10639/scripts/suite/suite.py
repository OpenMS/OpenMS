"""Frozen #10364 20-file suite for the total peptide order: arms old (develop merge-base) and new, each on the
normal database (PSMs) and on the entrapment database (FDP; Velos: natural entrapment on the normal run).
Arms interleaved per file. usage: suite.py <new build> [datasets]"""
import json, os, subprocess, sys, time
from pathlib import Path
S = Path('/home/user/bench/suite'); D = Path('/home/user/bench/frozen/z/prose-andes-reproduction/broad_benchmark/data')
DS = json.loads((S/'cfg/datasets.json').read_text())
ARMS = [('old', '/opt/openms-base'), ('new', sys.argv[1])]
for d in sys.argv[2:] or list(DS):
    mz = D/f'{d}.mzML'
    while not any(l.split() == [d, 'OK'] for f in S.glob('prep*.log') for l in open(f)):
        time.sleep(30)
    for db in ('normal', 'entrap'):
        fa = D/DS[d]['fasta'] if db == 'normal' else DS[d]['entrap'] and S/DS[d]['entrap']
        if not fa: continue
        for arm, b in ARMS:
            out = S/'runs'/d/f'{arm}_{db}'
            if (out/'metrics.json').exists(): continue
            print(d, arm, db, flush=True)
            subprocess.run(['python3', '-I', str(S/'scripts/run_arm.py'), '--bin', b, '--params', str(S/f'cfg/{d}.json'),
                            '--mzml', str(mz), '--fasta', str(fa), '--out', str(out)], check=True,
                           env=dict(os.environ, PERCOLATOR='/home/user/bench/tp/Linux/x86_64/Percolator/percolator'))
