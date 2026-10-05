import json,os,subprocess,time
from pathlib import Path
import argparse
parser=argparse.ArgumentParser(description="Compare Release identification I/O APIs on Linux")
parser.add_argument('binary',type=Path)
parser.add_argument('output',type=Path)
args=parser.parse_args()
binary=args.binary.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
env=os.environ.copy();env['OMP_NUM_THREADS']='1'
records=[]
def measure(label,args):
 print(label,flush=True)
 start=time.monotonic()
 with (out/(label+'.log')).open('w') as f:
  child=subprocess.Popen([str(binary)]+[str(x) for x in args],env=env,stdout=f,stderr=subprocess.STDOUT,text=True)
  pid,status,usage=os.wait4(child.pid,0);child.returncode=os.waitstatus_to_exitcode(status)
 text=(out/(label+'.log')).read_text()
 phases={};digests=[]
 for line in text.splitlines():
  if line.startswith('TIME\t'):
   _,phase,value=line.split('\t');phases[phase]=float(value)
  if line.startswith('DIGEST\t'):digests.append(line.split('\t')[1:])
 record={'label':label,'args':[str(x) for x in args],'code':child.returncode,'elapsed':time.monotonic()-start,'rss_kib':usage.ru_maxrss,'phases':phases,'digests':digests}
 records.append(record);(out/'measurements.json').write_text(json.dumps(records,indent=2))
 print(json.dumps(record),flush=True)
 if child.returncode: raise RuntimeError(label+' failed:\n'+text[-3000:])
 if not digests or any(d!=digests[0] for d in digests):raise RuntimeError('Conversion digest mismatch '+label)
 return digests[0]
for rows,runs in [(1000,2),(100000,1),(1000000,1),(1000000,1000)]:
 case=f'{rows}-{runs}';expected=None
 formats=['idxml','parquet','oms','native'];paths={f:out/(case+'-'+f) for f in formats}
 for f in formats:
  got=measure(case+'-'+f+'-write',['write',f,paths[f],rows,runs])
  if expected is None:expected=got
  if got!=expected:raise RuntimeError('Input mismatch '+f)
 for repetition in range(4):
  for f in formats[repetition:]+formats[:repetition]:
   got=measure(case+'-'+f+'-read-'+str(repetition),['read',f,paths[f]])
   if got!=expected:raise RuntimeError('Roundtrip mismatch '+f+': '+str(got)+' versus '+str(expected))
 sizes={f:sum(p.stat().st_size for p in paths[f].rglob('*') if p.is_file()) if paths[f].is_dir() else paths[f].stat().st_size for f in formats}
 (out/(case+'-sizes.json')).write_text(json.dumps(sizes,indent=2))
print('Comparison completed',flush=True)
