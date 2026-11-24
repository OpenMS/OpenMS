"""Read-only, standard-library identity/resource comparison. No scientific rescoring."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def artifacts(raw):
    stem = f'LFQ_Astral_DIA_15min_50ng_Human_{raw:02d}'
    return ['precursors.tsv', stem + '.precursors.tsv', stem + '.calibration.json']


def accounting(path):
    if not path.exists():
        return {}
    with path.open() as stream:
        return {row['JobID']: row for row in csv.DictReader(stream, delimiter='|')}


def resources(directory):
    path = directory / 'time.txt'
    text = path.read_text() if path.exists() else ''
    rss = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', text)
    wall = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)', text)
    seconds = None
    if wall:
        seconds = 0.0
        for part in wall[1].split(':'):
            seconds = seconds * 60 + float(part)
    return {'wall_seconds': seconds, 'max_rss_kib': int(rss[1]) if rss else None}


def compare(candidate, baseline, raw, job, accounts):
    exit_path = candidate / 'exit_code.txt'
    exit_code = exit_path.read_text().strip() if exit_path.exists() else 'MISSING'
    record = accounts.get(job, {})
    files = {}
    for name in artifacts(raw):
        a, b = candidate / name, baseline / name
        files[name] = ('MISSING' if not a.is_file() or not b.is_file() else
                       'SAME' if digest(a) == digest(b) else 'DIFFERENT')
    baseline_exit = baseline / 'exit_code.txt'
    baseline_ok = baseline_exit.exists() and baseline_exit.read_text().strip() == '0'
    success = (exit_code == '0' and baseline_ok and record.get('State') == 'COMPLETED'
               and record.get('ExitCode') == '0:0' and all(v == 'SAME' for v in files.values()))
    return {'status': 'PASS' if success else 'FAILED_OR_MISSING', 'exit_code': exit_code,
            'baseline_exit_ok': baseline_ok, 'job_id': job, 'accounting': record, 'files': files,
            'candidate': resources(candidate), 'baseline': resources(baseline)}


def evaluate(work, baseline, accounts, array_job):
    rows = []
    for i in range(15):
        raw, seed = i // 5 + 1, i % 5
        name = f'raw{raw:02d}_s{seed}'
        rows.append(dict(raw=raw, seed=seed, **compare(work / 'runs' / name,
                    baseline / 'runs' / name, raw, f'{array_job}_{i}', accounts)))
    return {'denominator': 15, 'passed': sum(r['status'] == 'PASS' for r in rows),
            'interpretation': 'Matched descriptive resources only. Hardware and cache variability apply. No FDR or quantification claim.',
            'runs': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('--bootstrap', action='store_true')
    args = parser.parse_args()
    accounts = accounting(args.work / ('bootstrap_accounting.tsv' if args.bootstrap else 'results/slurm_accounting.tsv'))
    if args.bootstrap:
        job = (args.work / 'bootstrap_job_id.txt').read_text().strip()
        result = compare(args.work / 'bootstrap', args.baseline / 'runs/raw01_s0', 1, job, accounts)
        print(json.dumps(result, indent=2))
        return 0 if result['status'] == 'PASS' else 1
    job = (args.work / 'primary_job_id.txt').read_text().strip()
    result = evaluate(args.work, args.baseline, accounts, job)
    print(json.dumps(result, indent=2))
    return 0 if result['passed'] == 15 else 1


if __name__ == '__main__':
    raise SystemExit(main())
