#!/usr/bin/env python3
"""Acquisition properties relevant to PR #10631, from scan_table.py output:
FAIMS CVs, how often the last MS1 in time is at another CV / is not the referenced survey scan,
isolation target != selected ion m/z, precursor outside its window, reporter scan windows above 126 Th."""
import sys, csv, collections
r = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
lowest_reporter = float(sys.argv[2]) if len(sys.argv) > 2 else 126.0
idx = {x['native_id']: i for i, x in enumerate(r)}
lv = collections.Counter(x['ms_level'] for x in r)
cvs = collections.Counter(x['faims_cv'] for x in r if x['ms_level'] == '1')
last = None; n2 = 0; refdiff = 0; lastcv = 0; noref = 0; offc = 0; outside = 0
for x in r:
    if x['ms_level'] == '1': last = x; continue
    if x['ms_level'] != '2': continue
    n2 += 1
    j = idx.get(x['spectrum_ref'])
    if j is None: noref += 1
    elif last is not None and r[j]['native_id'] != last['native_id']: refdiff += 1
    if last is not None and last['faims_cv'] != x['faims_cv']: lastcv += 1
    if x['iso_target'] and x['selected_mz']:
        t, s = float(x['iso_target']), float(x['selected_mz'])
        if abs(t - s) > 0.005: offc += 1
        if x['iso_lo'] and x['iso_hi'] and not (t - float(x['iso_lo']) <= s <= t + float(x['iso_hi'])): outside += 1
qlevel = '3' if lv.get('3') else '2'
q = [x for x in r if x['ms_level'] == qlevel]
above = sum(1 for x in q if x['sw_lo'] and float(x['sw_lo']) > lowest_reporter)
iso = collections.Counter(x['iso_lo'] for x in r if x['ms_level'] == '2').most_common(3)
print(f"{sys.argv[1]}\tlevels={dict(lv)}\tFAIMS_CVs(MS1)={dict(cvs)}\tMS2={n2}\tno_spectrumRef={noref}\tref!=lastMS1={refdiff}\tlastMS1_otherCV={lastcv}"
      f"\ttarget!=selected={offc}\tselected_outside_window={outside}\tiso_lower_offsets={iso}\treporter_scans(MS{qlevel})={len(q)}\tscan_window_start>{lowest_reporter}={above}")
