"""Run the frozen #10364 20-file suite, arms interleaved per file. Waits for inputs still being prepared."""
import json, subprocess, sys, time
from pathlib import Path
R = Path('/home/user/bench'); D = R/'frozen/broad_benchmark/data'
DS = json.loads((R/'cfg/datasets.json').read_text())
ON = 'Search:report:collapse_positional_isomers=true'
ARMS = [('base', '/opt/openms-base', False, []), ('off', '/opt/openms-pr', False, []), ('on', '/opt/openms-pr', False, [ON]),
        ('off_entrap', '/opt/openms-pr', True, []), ('on_entrap', '/opt/openms-pr', True, [ON])]
order = sys.argv[1:] or list(DS)
for d in order:
    mz = D/f'{d}.mzML'
    while not mz.exists() or 'OK' not in [l.split()[1] for l in open(R/'prep.log').readlines() + open(R/'prep_tims.log').readlines() if l.startswith(d + ' ') and len(l.split()) == 2] :
        time.sleep(30)
    for arm, b, entrap, extra in ARMS:
        out = R/'runs/suite'/d/arm
        if (out/'metrics.json').exists(): continue
        fa = DS[d]['fasta'].replace('_td', '_entrap_td') if entrap else DS[d]['fasta']
        cmd = ['python3', '-I', str(R/'scripts/run_arm.py'), '--bin', b, '--params', str(R/f'cfg/{d}.json'), '--mzml', str(mz),
               '--fasta', str(D/fa), '--out', str(out)] + sum([['--set', x] for x in extra], [])
        print(d, arm, flush=True)
        subprocess.run(cmd, check=True)
