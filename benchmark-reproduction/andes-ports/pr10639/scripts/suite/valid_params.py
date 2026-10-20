"""Write <bin>/valid_params.json: the Search: parameter names of <bin>/bin/ProSE -write_ini."""
import json, os, subprocess, sys, tempfile
import xml.etree.ElementTree as ET
b = sys.argv[1]; ini = tempfile.mktemp(suffix='.ini')
subprocess.run([f'{b}/bin/ProSE', '-write_ini', ini], check=True, env=dict(os.environ, OPENMS_DATA_PATH=f'{b}/share/OpenMS', LD_LIBRARY_PATH=f'{b}/lib'))
names = []
def walk(n, path):
    for c in n:
        if c.tag == 'NODE': walk(c, path + c.get('name') + ':')
        elif c.tag in ('ITEM', 'ITEMLIST') and path.startswith('ProSE:1:Search:'): names.append(path[len('ProSE:1:Search:'):] + c.get('name'))
walk(ET.parse(ini).getroot(), ''); os.unlink(ini)
json.dump(names, open(f'{b}/valid_params.json', 'w')); print(b, len(names))
