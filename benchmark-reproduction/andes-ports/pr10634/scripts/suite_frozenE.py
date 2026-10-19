"""Frozen #10364 paired-mass entrapment (r = 1; sage-comparison-20261002/make_entrapment.py, sha256 fb77e3df.../76d25071...)
for HF-X, Astral and TMTpro, PR build option off vs on."""
import json, subprocess
from pathlib import Path
R = Path('/home/user/bench'); D = R/'frozen/broad_benchmark/data'
DS = json.loads((R/'cfg/datasets.json').read_text())
ON = 'Search:report:collapse_positional_isomers=true'
for d in ['hfx_A2', 'hfx_B1', 'hfx_B3', 'astral_A2', 'astral_B1', 'astral_B3', 'eclipse_tmtpro_10855', 'eclipse_tmtpro_10858', 'eclipse_tmtpro_10863']:
    for arm, extra in (('off_frozenE', []), ('on_frozenE', ['--set', ON])):
        out = R/'runs/suite'/d/arm
        if (out/'metrics.json').exists(): continue
        fa = DS[d]['fasta'].replace('_td', '_frozenE_td')
        subprocess.run(['python3', '-I', str(R/'scripts/run_arm.py'), '--bin', '/opt/openms-pr', '--params', str(R/f'cfg/{d}.json'),
                        '--mzml', str(D/f'{d}.mzML'), '--fasta', str(D/fa), '--out', str(out)] + extra, check=True)
