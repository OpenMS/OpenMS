#!/usr/bin/env python3
"""Validate and summarize raw layout measurements without inventing workload weights."""
# $Maintainer: Justin Sing $
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import statistics


def require(condition, message):
    if not condition:
        raise ValueError(message)


def summarize(directory):
    manifest = json.loads((directory / 'manifest.json').read_text())
    require(manifest.get('completed_cases') == len(manifest['schedule']), 'Incomplete measurement matrix')
    rows = [json.loads(line) for line in (directory / 'runs.jsonl').read_text().splitlines()]
    require(len(rows) == manifest['completed_cases'], 'Missing or extra result records')
    expected = {(job['fixture'], job['min_points'], job['variant'], job['trial']) for job in manifest['schedule']}
    observed = {(row['fixture'], row['min_points'], row['variant'], row['trial']) for row in rows}
    require(len(expected) == len(rows) and expected == observed, 'Duplicate, missing or unexpected cases')
    groups = defaultdict(list)
    checksums = {}
    common = {}
    workloads = {}
    for row in rows:
        require(row['returncode'] == 0, 'Failed benchmark process')
        result = row['result']
        require(result['schema_version'] == 1 and result['variant'] == row['variant'], 'Incorrect output schema/variant')
        require(result['rows'] == result['source_rows'] * result['replication_factor'], 'Row replication mismatch')
        require(result['scans'] == result['source_scans'] * result['replication_factor'], 'Scan replication mismatch')
        require(result['rows'] >= row['min_points'], 'Requested scale not met')
        require(result['allocator_delta_bytes'] > 0, 'Nonpositive retained allocation delta')
        key = (row['fixture'], row['min_points'])
        geometry = tuple(result[name] for name in ('rows', 'scans', 'source_rows', 'source_scans', 'targets'))
        require(common.setdefault(key, geometry) == geometry, 'Variants measured different inputs')
        require(workloads.setdefault(key, set(result['timings'])) == set(result['timings']),
                'Variants measured different workloads')
        for name, timing in result['timings'].items():
            require(len(timing['samples_ms']) == manifest['repeats'], 'Incorrect timing repeat count')
            require(all(math.isfinite(value) and value >= 0 for value in timing['samples_ms']), 'Invalid timing sample')
            require(math.isfinite(timing['checksum']), 'Nonfinite timing checksum')
            checksum_key = key + (name,)
            reference = checksums.setdefault(checksum_key, timing['checksum'])
            require(math.isclose(reference, timing['checksum'], rel_tol=1e-10, abs_tol=1e-8),
                    f'Cross-variant/trial checksum mismatch: {checksum_key}')
        groups[key + (row['variant'],)].append(result)
    summary = []
    for (fixture, scale, variant), results in sorted(groups.items()):
        require(len(results) == manifest['trials'], 'Unequal outer trial counts')
        first = results[0]
        require(all(set(result['timings']) == set(first['timings']) for result in results), 'Workload set changed')
        item = dict(fixture=fixture, min_points=scale, variant=variant, trials=len(results),
                    rows=first['rows'], scans=first['scans'], source_rows=first['source_rows'],
                    source_scans=first['source_scans'], replication_factor=first['replication_factor'])
        for name in ('memory_snapshot_stage', 'build_validation', 'public_kernel_adapter', 'memory_scope',
                     'public_reference', 'sizeof_raw_point', 'sizeof_coordinates', 'sizeof_hot_point',
                     'sizeof_scan_metadata', 'sparse_rows', 'fullrow_rows'):
            if name in first:
                item[name] = first[name]
        for name in ('allocator_delta_bytes', 'baseline_rss_bytes', 'rss_bytes', 'peak_rss_bytes',
                     'construction_rss_bytes', 'construction_peak_rss_bytes', 'final_peak_rss_bytes', 'build_ms'):
            if name in first:
                values = [result[name] for result in results]
                item[name] = dict(median=statistics.median(values), min=min(values), max=max(values))
        item['bytes_per_point'] = item['allocator_delta_bytes']['median'] / item['rows']
        item['timings'] = {}
        for name in sorted(first['timings']):
            outer = [result['timings'][name]['median_ms'] for result in results]
            samples = [value for result in results for value in result['timings'][name]['samples_ms']]
            item['timings'][name] = dict(median_ms=statistics.median(outer), outer_min_ms=min(outer),
                                          outer_max_ms=max(outer), sample_min_ms=min(samples),
                                          sample_max_ms=max(samples), samples_ms=samples)
        summary.append(item)
    return dict(schema_version=1, completed_cases=len(rows), source_commit=manifest['source_commit'],
                source_sha256=manifest['source_sha256'], fixture_sha256=manifest['fixture_sha256'],
                binary_sha256=manifest['binary_sha256'], started_utc=manifest['started_utc'],
                finished_utc=manifest['finished_utc'], cpu=manifest['cpu'], governor=manifest['governor'],
                last_level_cache=manifest['last_level_cache'], seed=manifest['seed'],
                repeats=manifest['repeats'], trials=manifest['trials'], cases=summary,
                interpretation='Median of fresh-process medians. Build validation differs across candidates. '
                               'Live allocator delta is not payload-only or process RSS. '
                               'Scaled data repeats source acquisitions. No whole-engine/global-optimal claim.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.directory)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(f'Validated {result["completed_cases"]} isolated runs and {len(result["cases"])} summarized cases')
    for case in result['cases']:
        times = ' '.join(f'{name}={data["median_ms"]:.4f}ms' for name, data in case['timings'].items())
        print(f'{case["fixture"]:5s} {case["rows"]:8d} {case["variant"]:12s} '
              f'{case["bytes_per_point"]:7.2f} B/point {times}')


if __name__ == '__main__':
    main()
