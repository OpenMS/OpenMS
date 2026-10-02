"""How much of PeptDeepRescoring's calibration set (top half of all hits by search score) is rank-1 target."""
import gzip, shutil, sys, statistics, tempfile
import pyopenms as p
for d in sys.argv[1:]:
    src = f"results/{d}/pd_on/search.idXML.gz"
    with tempfile.NamedTemporaryFile(suffix=".idXML") as tmp:
        with gzip.open(src, "rb") as f: shutil.copyfileobj(f, tmp)
        tmp.flush()
        prots, peps = [], p.PeptideIdentificationList()
        p.IdXMLFile().load(tmp.name, prots, peps)
    rows = []
    for pi in peps:
        for r, h in enumerate(pi.getHits()):
            rows.append((h.getScore(), r, h.getMetaValue("target_decoy"), float(h.getMetaValue("rt_abs_error")), pi.getRT()))
    rows.sort(key=lambda x: -x[0])
    conf = rows[: len(rows) // 2]
    n = len(conf)
    r1t = [x for x in conf if x[1] == 0 and x[2] == "target"]
    print(d, "hits", len(rows), "confident", n, "rank1 target", len(r1t), f"({100*len(r1t)/n:.1f}%)",
          "decoy", sum(x[2] == "decoy" for x in conf), "rank>1", sum(x[1] > 0 for x in conf))
    for lab, sel in (("rank1 target", [x for x in rows if x[1] == 0 and x[2] == "target"]),
                     ("rank1 decoy", [x for x in rows if x[1] == 0 and x[2] == "decoy"]),
                     ("rank>1", [x for x in rows if x[1] > 0])):
        print("   ", lab, len(sel), "median rt_abs_error", round(statistics.median(x[3] for x in sel), 1))
    top = [x for x in rows if x[1] == 0 and x[2] == "target"][:1000]
    print("    top-1000 rank1 targets median rt_abs_error", round(statistics.median(x[3] for x in top), 1), "RT range", round(min(x[4] for x in rows)), round(max(x[4] for x in rows)))
