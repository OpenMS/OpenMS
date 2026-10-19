"""Stream an mzML from stdin and keep N MS2 spectra with a selected-ion m/z, chosen like the #10364 prepare.py:
the N smallest SHA256(dataset ID + NUL + native ID), written in their original (RT) order, non-indexed mzML.
usage: curl ... | subsample_mzml.py <dataset id> <out.mzML> [N=8000]"""
import hashlib, heapq, re, sys
ds, out = sys.argv[1], sys.argv[2]; N = int(sys.argv[3]) if len(sys.argv) > 3 else 8000
buf = b""; header = None; heap = []; idx = 0; total = 0
src = sys.stdin.buffer
while True:
    chunk = src.read(1 << 22)
    if chunk: buf += chunk
    if header is None:
        i = buf.find(b"<spectrumList")
        if i < 0:
            if not chunk: break
            continue
        header = buf[:i]; j = buf.find(b">", i); buf = buf[j + 1:]
    while True:
        s = buf.find(b"<spectrum ")
        if s < 0: break
        e = buf.find(b"</spectrum>", s)
        if e < 0: break
        el = buf[s:e + len(b"</spectrum>")]; buf = buf[e + len(b"</spectrum>"):]
        idx += 1
        if b'accession="MS:1000511" name="ms level" value="2"' not in el or b'MS:1000744' not in el: continue
        total += 1
        nid = re.search(rb'id="([^"]*)"', el).group(1)
        h = hashlib.sha256(ds.encode() + b"\0" + nid).hexdigest()
        item = (tuple(-ord(c) for c in h[:16]), h, idx, el)  # max-heap on hash via negation
        if len(heap) < N: heapq.heappush(heap, item)
        elif item > heap[0]: heapq.heapreplace(heap, item)
    if not chunk: break
sel = sorted(heap, key=lambda x: x[2])
with open(out, "wb") as f:
    hdr = re.sub(rb"<indexedmzML[^>]*>", b"", header)
    f.write(hdr + b'<spectrumList count="%d" defaultDataProcessingRef="' % len(sel))
    m = re.search(rb'<spectrumList[^>]*defaultDataProcessingRef="([^"]*)"', buf + header)
    dp = re.search(rb'<dataProcessing id="([^"]*)"', header)
    f.write((dp.group(1) if dp else b"dp") + b'">\n')
    for n, (_, _, _, el) in enumerate(sel):
        f.write(re.sub(rb'index="\d+"', b'index="%d"' % n, el, count=1) + b"\n")
    f.write(b"</spectrumList>\n</run>\n</mzML>\n")
print("MS2 with precursor:", total, "selected:", len(sel), file=sys.stderr)
