"""Per file and seed: change in accepted PSMs and entrapment hits of an arm against develop+fix
(hds1_E / tds1_E), and a sign test over the (file, seed) pairs whose entrapment count changed."""
import math
from pathlib import Path
R = Path(__file__).resolve().parent.parent / "results"
FILES = {"Astral HCD": ["astral_A2", "astral_B1", "astral_B3"],
         "Exploris 480 TMTpro": ["eclipse_tmtpro_10855", "eclipse_tmtpro_10858", "eclipse_tmtpro_10863"],
         "HF-X HCD": ["hfx_A2", "hfx_B1", "hfx_B3"]}


def acc(d, a, s):
    out = {}
    for l in (R / d / a / f"s{s}" / "target.tsv").read_text().splitlines()[1:]:
        f = l.split("\t")
        if float(f[2]) <= 0.01:
            out[f[0].rsplit("_", 2)[-2]] = all(x.startswith("ENTRAP_") for x in f[5:] if x)
    return out


def sign_p(pos, neg):
    n, k = pos + neg, min(pos, neg)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n) if n else 1.0


for arm in ("d_priors_E", "d_79_E"):
    for g, files in FILES.items():
        cells, pos, neg = [], 0, 0
        for d in files:
            base = "tds1_E" if d.startswith("eclipse") else "hds1_E"
            for s in (1, 42, 137):
                a, b = acc(d, arm, s), acc(d, base, s)
                de = sum(a.values()) - sum(b.values())
                pos, neg = pos + (de > 0), neg + (de < 0)
                cells.append(f"{len(a) - len(b):+d}/{de:+d}")
        print(f"{arm:11s} {g:20s} entrapment up in {pos}, down in {neg} of 9 (sign test p = {sign_p(pos, neg):.3f}); "
              f"PSMs/entrapment per file x seed: {' '.join(cells)}")
