"""List accepted PSMs where protein-based (EntrapmentEvaluator rule) and sequence-based (frozen) Velos entrapment disagree."""
import csv, sys, re
sys.path.insert(0, sys.argv[1]); import frozen_eval as fe
entries = fe.read_fasta(sys.argv[2]); db, _ = fe.classes(entries)
raw = {"sample": [], "ups1": [], "human": []}  # without I->L
for h, s in entries:
    acc = h.split("|")[1] if "|" in h else h.split()[0]
    k = "sample" if ("Cont_" in h or "_YEAST" in h.split()[0]) else "ups1" if acc in fe.UPS1 else "human" if "_HUMAN" in h.split()[0] else "sample"
    raw[k].append(s)
raw = {k: "\n".join(v) for k, v in raw.items()}
for path in sys.argv[3:]:
    for row in open(path).read().splitlines()[1:]:
        f = row.split("\t")
        if float(f[2]) > 0.01: continue
        prot_e = all(p.startswith("ENTRAP_") for p in f[5:] if p)
        pep = fe.stripped(f[4]); cat = "sample" if pep in db["sample"] else "ups1" if pep in db["ups1"] else "entrapment" if pep in db["human"] else "unknown"
        if prot_e != (cat == "entrapment"):
            plain = re.sub(r"[^A-Z]", "", re.sub(r"\[[^\]]*\]", "", f[4]).split(".")[1])
            print(path.split("/")[-2], f[4], cat, "exact-in-sample" if plain in raw["sample"] else "-", "exact-in-ups1" if plain in raw["ups1"] else "-", ";".join(f[5:])[:120])
