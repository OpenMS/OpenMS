"""Median wall time / peak RSS per experiment and arm (all reps), and whether idXML/PIN equal the reference arm.
usage: summarize.py <results/*.tsv>  -> markdown table on stdout"""
import csv, statistics, sys
from collections import defaultdict
print('| experiment | arm | threads | runs | wall s (median) | peak RSS GiB (median) | idXML = ref | PIN = ref | status |')
print('|---|---|---:|---:|---:|---:|---|---|---|')
for path in sys.argv[1:]:
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    by = defaultdict(list)
    for r in rows: by[r['arm']].append(r)
    ref_arm = 's1' if 's1' in by else next(iter(by))
    ref = next((r for r in by[ref_arm] if r['status'] == 'ok'), None)
    for arm, rs in by.items():
        ok = [r for r in rs if r['status'] == 'ok']
        wall = statistics.median(float(r['wall_s']) for r in ok) if ok else float('nan')
        rss = statistics.median(int(r['maxrss_kb']) for r in ok) / 2**20 if ok else max(int(r['maxrss_kb'] or 0) for r in rs) / 2**20
        same = lambda k: ('yes' if all(r[k] == ref[k] for r in ok) else 'NO') if ok and ref else '-'
        st = ','.join(sorted({r['status'] for r in rs}))
        print(f"| {rows[0]['exp']} | {arm} | {rs[0]['threads']} | {len(rs)} | {wall:.1f} | {rss:.2f} | {same('idxml_sha256')} | {same('pin_sha256')} | {st} |")
