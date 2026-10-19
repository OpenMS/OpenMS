# id_bruker (PXD058436, MHCquant2 timsTOF Pro, class I): 8,000 MS2 per file with the frozen #10364 prepare.py sampling
set -e
I=/home/user/bench/idb; cd $I
[ -f data/human_sp.fasta ] || curl -sSfL -o data/human_sp.fasta https://archive.openms.de/openms/testfiles/human_sp.fasta
for d in idb_spleen idb_corpus idb_liver; do
  [ -f data/$d.mzML ] && continue
  url=$(python3 -c "import json;print(next(r['url'] for r in json.load(open('manifest.json')) if r['id']=='$d'))")
  size=$(python3 -c "import json;print(next(r['bytes'] for r in json.load(open('manifest.json')) if r['id']=='$d'))")
  for i in $(seq 1 30); do
    [ -f raw/$d.mzML ] && [ $(stat -c %s raw/$d.mzML) -ge $size ] && break
    curl -sSL -C - -o raw/$d.mzML "$url" || true; sleep 3
  done
  [ $(stat -c %s raw/$d.mzML) -eq $size ] || { echo "size mismatch $d"; exit 1; }
  ls -l raw/$d.mzML; sha256sum raw/$d.mzML
  python3 prepare.py $d
  sha256sum data/$d.mzML
  rm -f raw/$d.mzML
done
echo IDB_INPUTS_DONE
