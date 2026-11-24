#!/usr/bin/env python3
"""Read-only monitor for the three frozen Codon job IDs. Never submits or cancels."""
import datetime
import json
from pathlib import Path
import shlex
import subprocess
import time

WORK = '/nfs/research/juan/singj/opendia/astral-parquet-20261124'
IDS = '42109936,42109937,42109938'
EVIDENCE = Path('/home/singjc/.jcode/scratch/compact-20261124/codon-monitor')
EVIDENCE.mkdir(exist_ok=True)
REMOTE = f'''export PATH=/ebi/slurm/codon/bin:$PATH
printf '__QUEUE__\\n'
squeue -h -j {IDS} -o '%i|%T|%M|%R'
printf '__ACCOUNTING__\\n'
sacct -X -n -P -j {IDS} --format=JobID,State,ExitCode,Elapsed,NodeList
printf '__BOOTSTRAP_LOG__\\n'
if test -f {WORK}/bootstrap/run.log; then tail -4 {WORK}/bootstrap/run.log; else tail -4 {WORK}/logs/bootstrap-* 2>/dev/null || true; fi
printf '__EVALUATION__\\n'
if test -f {WORK}/results/evaluation_exit_codes.txt; then cat {WORK}/results/evaluation_exit_codes.txt; fi
'''
SSH = ['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20', '-o', 'ControlMaster=no',
       '-o', 'ControlPath=none', '-o', 'ServerAliveInterval=30', '-o', 'ServerAliveCountMax=2',
       'codon', 'become pst_prd bash -lc ' + shlex.quote(REMOTE)]
terminal = {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL',
            'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED'}
for attempt in range(288):
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    try:
        result = subprocess.run(SSH, text=True, capture_output=True, timeout=120)
        output = result.stdout
        (EVIDENCE / f'{attempt:03d}.txt').write_text(now + '\n' + output + '\nSTDERR:\n' + result.stderr)
        if result.returncode or '__ACCOUNTING__' not in output:
            message = f'Codon read-only monitor transport/query issue {result.returncode}, no workload action'
        else:
            queue = output.split('__QUEUE__\n', 1)[1].split('__ACCOUNTING__\n', 1)[0].strip()
            accounting = output.split('__ACCOUNTING__\n', 1)[1].split('__BOOTSTRAP_LOG__\n', 1)[0].strip()
            print(now + '\n' + output, flush=True)
            message = 'Codon frozen experiment active: ' + ('; '.join(queue.splitlines()) or 'queue empty, awaiting accounting')
            records = [[field.strip() for field in line.split('|')] for line in accounting.splitlines()]
            evaluation = next((row for row in records if row[0] == '42109938'), None)
            bootstrap = next((row for row in records if row[0] == '42109936'), None)
            if not queue and evaluation and bootstrap and all(row[1].split()[0] in terminal for row in [evaluation, bootstrap]):
                print('JCODE_PROGRESS ' + json.dumps({'message': 'Frozen Codon jobs terminal, collect fixed-denominator evaluation'}), flush=True)
                raise SystemExit(0)
    except (subprocess.TimeoutExpired, OSError) as error:
        message = f'Codon monitor read-only query failed: {error}'
        (EVIDENCE / f'{attempt:03d}.error.txt').write_text(now + '\n' + str(error))
    print('JCODE_PROGRESS ' + json.dumps({'message': message}), flush=True)
    time.sleep(300)
raise SystemExit('Read-only monitor reached 24-hour active polling limit, inspect existing jobs without resubmission')
