"""Write a ProSE INI with the frozen #10364 native_dedup parameters of one dataset (jobs.json), plus overrides.
usage: gen_ini.py <template.ini> <jobs.json> <dataset> <out.ini> [key=value ...]   (key without 'Search:'; lists with ',')"""
import json, sys
import xml.etree.ElementTree as ET
tmpl, jobs, dataset, out = sys.argv[1:5]
LISTS = ('modifications:fixed', 'modifications:variable')
job = next(j for j in json.load(open(jobs)) if j['id'] == f'native_dedup:{dataset}:default')
p = dict(job['parameters'])
p.pop('snes_enabled', None); p.pop('precursor:isotopes', None)   # not ProSE TOPP parameters
for kv in sys.argv[5:]:
    k, v = kv.split('=', 1)
    p[k] = v.split(',') if (k in LISTS and v) else ([] if k in LISTS else v)
t = ET.parse(tmpl); root = t.getroot()
items = {}
def walk(n, path):
    for c in n:
        if c.tag == 'NODE': walk(c, path + c.get('name') + ':')
        elif c.tag in ('ITEM', 'ITEMLIST'): items[path + c.get('name')] = c
walk(root, '')
for k, v in p.items():
    key = 'ProSE:1:Search:' + k
    if key not in items: sys.exit('unknown parameter ' + k)
    e = items[key]
    if e.tag == 'ITEMLIST':
        for c in list(e): e.remove(c)
        for x in (v if isinstance(v, list) else [v]): ET.SubElement(e, 'LISTITEM', value=str(x))
    else:
        e.set('value', str(v).lower() if isinstance(v, bool) else str(v))
t.write(out)
