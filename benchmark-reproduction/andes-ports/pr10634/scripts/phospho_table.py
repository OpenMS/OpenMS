"""Phospho (PXD007653) table from runs/phospho/{base,off,on}/metrics.json."""
import json, sys
from pathlib import Path
R = Path(sys.argv[1])
rows = ['arm\tseed\tpsms\tpeptides\tpeptidoforms\tphospho_psms\tentrapment_psms\tfdp_pct\tpin_rows_target\tpin_rows_decoy\t'
        'isomer_sets_target\tisomer_sets_decoy\tisomer_extra_rows_target\tisomer_extra_rows_decoy\tpin_sha256']
for a in ('base', 'off', 'on'):
    m = json.loads((R/a/'metrics.json').read_text()); p = m['pin']
    for s, v in m['seeds'].items():
        rows.append('\t'.join(map(str, [a, s, v['psms'], v['peptides'], v['peptidoforms'], v['phospho_psms'], v['entrapment_psms'],
                    f"{100 * v['fdp']:.2f}", p['rows_target'], p['rows_decoy'], p['isomer_sets_target'], p['isomer_sets_decoy'],
                    p['isomer_extra_rows_target'], p['isomer_extra_rows_decoy'], m['pin_sha256'][:16]])))
print('\n'.join(rows))
