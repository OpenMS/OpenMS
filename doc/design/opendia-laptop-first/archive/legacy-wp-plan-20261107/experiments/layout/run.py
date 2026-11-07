#!/usr/bin/env python3
"""Run isolated, randomized layout comparisons. No third-party Python packages."""
# $Maintainer: Justin Sing $
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

VARIANTS = ('rawtile', 'aos', 'soa', 'factored-aos', 'factored-soa')
FIXTURES = {
    'nonim': 'topp/OpenSwathWorkflow_1_input.mzML',
    'im': 'topp/OpenSwathWorkflow_23_input.mzML',
    'faims': 'class_tests/openms/data/IM_FAIMS_test.mzML',
}


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def parse_result(stdout):
    candidates = []
    for line in stdout.splitlines():
        line = line.strip().removeprefix('LAYOUT_JSON ').removeprefix('RESULT ').strip()
        if line.startswith('{'):
            try:
                candidates.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    if len(candidates) != 1 or not isinstance(candidates[0], dict):
        raise ValueError('Expected exactly one JSON object from benchmark')
    return candidates[0]


def jobs(sizes, trials, large_im, seed):
    result = [dict(fixture=fixture, min_points=size, variant=variant, trial=trial)
              for trial, fixture, size, variant in itertools.product(
                  range(trials), FIXTURES, sizes, VARIANTS)]
    if large_im and large_im not in sizes:
        result.extend(dict(fixture='im', min_points=large_im, variant=variant, trial=trial)
                      for trial, variant in itertools.product(range(trials), VARIANTS))
    random.Random(seed).shuffle(result)
    return result


def self_test():
    assert parse_result('loader message\nRESULT {"ok": true}\n') == {'ok': True}
    assert parse_result('  LAYOUT_JSON {"schema_version": 1}\n') == {'schema_version': 1}
    try:
        parse_result('{"one": 1}\n{"two": 2}')
    except ValueError:
        pass
    else:
        raise RuntimeError('Ambiguous benchmark results accepted')
    first = jobs([0, 100], 3, 1000, 17)
    assert len(first) == 105 and first == jobs([0, 100], 3, 1000, 17)
    assert len({tuple(sorted(job.items())) for job in first}) == 105
    assert first != jobs([0, 100], 3, 1000, 18)
    print('Runner self-check passed: parsing, ambiguity, balanced unique jobs and deterministic randomization')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='New evidence directory, never overwritten')
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--trials', type=int, default=3)
    parser.add_argument('--repeats', type=int, default=7)
    parser.add_argument('--sizes', default='0,1048576')
    parser.add_argument('--large-im', type=int, default=4194304)
    parser.add_argument('--seed', type=int, default=20261025)
    args = parser.parse_args()
    sizes = [int(size) for size in args.sizes.split(',')]
    if args.trials < 1 or args.repeats < 1 or not sizes or min(sizes) < 0 or args.large_im < 0:
        parser.error('Trials/repeats must be positive and point scales nonnegative')
    if len(set(sizes)) != len(sizes):
        parser.error('Point scales must not contain duplicates')
    if args.cpu not in os.sched_getaffinity(0):
        parser.error('Requested CPU is outside allowed affinity')
    args.binary = args.binary.resolve(strict=True)
    args.source = args.source.resolve(strict=True)
    fixtures = {name: args.source / 'src/tests' / relative for name, relative in FIXTURES.items()}
    for path in fixtures.values():
        if not path.is_file():
            parser.error(f'Missing fixture {path}')
    schedule = jobs(sizes, args.trials, args.large_im, args.seed)
    args.output.mkdir(parents=True, exist_ok=False)
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    sources = {
        'records.cpp': args.source / 'src/topp/OpenDIA/source/OpenDIARecords.cpp',
        'records.h': args.source / 'src/topp/OpenDIA/include/OpenDIARecords.h',
        'extractor.cpp': args.source / 'src/openms/source/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.cpp',
        'benchmark.cpp': Path(__file__).with_name('layout_benchmark.cpp'),
        'CMakeLists.txt': Path(__file__).with_name('CMakeLists.txt'),
        'run.py': Path(__file__),
        'summarize.py': Path(__file__).with_name('summarize.py'),
    }
    manifest = {
        'started_utc': datetime.now(timezone.utc).isoformat(), 'binary_sha256': digest(args.binary),
        'source_commit': subprocess.check_output(['git', '-C', str(args.source), 'rev-parse', 'HEAD'], text=True).strip(),
        'source_dirty': subprocess.check_output(['git', '-C', str(args.source), 'status', '--porcelain'], text=True),
        'source_sha256': {name: digest(path) for name, path in sources.items()},
        'fixture_sha256': {name: digest(path) for name, path in fixtures.items()},
        'platform': platform.platform(), 'python': platform.python_version(), 'cpu': args.cpu,
        'lscpu': subprocess.check_output(['lscpu'], text=True),
        'threads': {key: environment[key] for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
        'seed': args.seed, 'trials': args.trials, 'repeats': args.repeats, 'schedule': schedule,
        'memory_note': 'Live allocator delta is distinct from RSS. Same source retained in every child.',
        'scale_note': 'Nonzero min_points repeats whole real fixture acquisitions, not new scientific data.',
        'host_note': 'One pinned CPU, not an exclusive host. No power/governor changes.',
    }
    cpu_directory = Path(f'/sys/devices/system/cpu/cpu{args.cpu}')
    for label, relative in (
        ('governor', 'cpufreq/scaling_governor'), ('cpu_max_khz', 'cpufreq/cpuinfo_max_freq'),
        ('siblings', 'topology/thread_siblings_list'), ('last_level_cache', 'cache/index3/size'),
    ):
        path = cpu_directory / relative
        manifest[label] = path.read_text().strip() if path.is_file() else None
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    with (args.output / 'runs.jsonl').open('x') as records:
        for index, job in enumerate(schedule):
            command = ['taskset', '-c', str(args.cpu), str(args.binary), str(fixtures[job['fixture']]),
                       job['variant'], str(job['min_points']), str(args.repeats)]
            started = time.monotonic()
            completed = subprocess.run(command, capture_output=True, text=True, errors='replace',
                                       env=environment, timeout=180)
            prefix = args.output / f'{index:03d}-{job["fixture"]}-{job["min_points"]}-{job["variant"]}'
            prefix.with_suffix('.stdout').write_text(completed.stdout)
            prefix.with_suffix('.stderr').write_text(completed.stderr)
            record = dict(job, command=command, returncode=completed.returncode,
                          process_seconds=time.monotonic() - started)
            if completed.returncode:
                records.write(json.dumps(record) + '\n')
                records.flush()
                raise RuntimeError(f'Benchmark failed, inspect {prefix}.stderr')
            record['result'] = parse_result(completed.stdout)
            records.write(json.dumps(record) + '\n')
            records.flush()
            print('JCODE_PROGRESS ' + json.dumps(dict(current=index + 1, total=len(schedule), unit='cases',
                                                       message=f'{job["fixture"]} {job["variant"]} {job["min_points"]}')), flush=True)
    # Reject changing source/binary/fixtures rather than mixing measurements across revisions.
    if digest(args.binary) != manifest['binary_sha256']:
        raise RuntimeError('Binary changed during measurements')
    for name, path in sources.items():
        if digest(path) != manifest['source_sha256'][name]:
            raise RuntimeError(f'Source changed during measurements: {name}')
    for name, path in fixtures.items():
        if digest(path) != manifest['fixture_sha256'][name]:
            raise RuntimeError(f'Fixture changed during measurements: {name}')
    manifest['finished_utc'] = datetime.now(timezone.utc).isoformat()
    manifest['completed_cases'] = len(schedule)
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Completed {len(schedule)} fresh-process comparisons. Raw evidence: {args.output}', flush=True)


if __name__ == '__main__':
    if sys.argv[1:] == ['--self-test']:
        self_test()
    else:
        main()
