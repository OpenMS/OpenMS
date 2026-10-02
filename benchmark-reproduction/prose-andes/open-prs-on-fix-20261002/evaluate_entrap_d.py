"""Entrapment FDP of the 2026-10-02 develop+fix arms (shuffled paired-mass entrapment, r = 1).

A PSM is an entrapment hit when every protein it maps to is ENTRAP_*. Per file: mean over seeds 1/42/137
of accepted target PSMs (q <= 0.01) and of entrapment hits. Per group: N and N_E are those seed means
summed over the group's files; combined FDP = 2 N_E / N. Baseline: develop with the deisotoping fix
(hds1_E / tds1_E; d_base_E on hfx_A2 checks that it equals develop 3a47278).
"""
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GROUPS = {"HF-X HCD": ["hfx_A2", "hfx_B1", "hfx_B3"], "Astral HCD": ["astral_A2", "astral_B1", "astral_B3"],
          "Exploris 480 TMTpro": ["eclipse_tmtpro_10855", "eclipse_tmtpro_10858", "eclipse_tmtpro_10863"]}
ARMS = {"develop + fix": lambda d: "tds1_E" if d.startswith("eclipse") else "hds1_E",
        "#10378 priors": lambda d: "d_priors_E", "#10379 default": lambda d: "d_79_E"}


def counts(d, arm):
    n, e = [], []
    for s in (1, 42, 137):
        p = ROOT / "results" / d / arm / f"s{s}" / "target.tsv"
        if not p.exists():
            return None
        acc = [l.split("\t") for l in p.read_text().splitlines()[1:]]
        acc = [f for f in acc if float(f[2]) <= 0.01]
        n.append(len(acc))
        e.append(sum(all(x.startswith("ENTRAP_") for x in f[5:] if x) for f in acc))
    return statistics.mean(n), statistics.mean(e)


def main():
    lines = ["group\tarm\tpsms\tentrapment\tcombined_fdp_pct"]
    for g, files in GROUPS.items():
        for name, arm in ARMS.items():
            c = [counts(d, arm(d)) for d in files]
            if any(x is None for x in c):
                lines.append(f"{g}\t{name}\tmissing\t\t")
                continue
            n, e = sum(x[0] for x in c), sum(x[1] for x in c)
            lines.append(f"{g}\t{name}\t{n:.1f}\t{e:.1f}\t{200 * e / n:.2f}")
    chk = (counts("hfx_A2", "d_base_E"), counts("hfx_A2", "hds1_E"))
    lines.append(f"# check hfx_A2: d_base_E {chk[0]} vs hds1_E {chk[1]}")
    text = "\n".join(lines) + "\n"
    (ROOT / "entrapment_d.tsv").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
