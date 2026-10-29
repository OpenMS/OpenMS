"""A small sequence-tag search used to test hypotheses before the real search.

De novo sequence tags (>= 3 residues) are read from the MS2 spectra and
matched, together with the precursor mass and the flanking fragment mass,
against in-silico digests of candidate databases. Target and reversed-decoy
digests are searched the same way, so every count comes with its own
background. The same machinery answers several Auto Mode questions:

* organism: which model-organism proteome explains the most spectra
  (net target - decoy matches; spectra that only one organism explains);
* Cys alkylation: which Cys mass explains the most Cys-containing peptides;
* protease: which cleavage rule explains the most spectra;
* mass accuracy: the precursor errors of the matched spectra give the
  systematic offset and the spread, i.e. a calibration search.

This is a prototype in pure Python/numpy. The production version belongs in
ProSE's own fragment index (a subset search per hypothesis), see the report.
"""

from __future__ import annotations

import gzip
import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

PROTON = 1.007276466812
H2O = 18.010564684
ISO = 1.0033548378

AA = {
    "G": 57.021464, "A": 71.037114, "S": 87.032028, "P": 97.052764, "V": 99.068414,
    "T": 101.047679, "C": 103.009185, "L": 113.084064, "I": 113.084064, "N": 114.042927,
    "D": 115.026943, "Q": 128.058578, "K": 128.094963, "E": 129.042593, "M": 131.040485,
    "H": 137.058912, "F": 147.068414, "R": 156.101111, "Y": 163.063329, "W": 186.079313,
}

CYS_MODS = {
    "Carbamidomethyl (C)": 57.021464,
    "none (free Cys)": 0.0,
    "Propionamide (C)": 71.037114,
    "Nethylmaleimide (C)": 125.047679,
    "Methylthio (C)": 45.987721,
    "Carboxymethyl (C)": 58.005479,
}

LABEL_MASS = {"TMT": 229.162932, "TMTpro": 304.207146, "iTRAQ4plex": 144.102063,
              "iTRAQ8plex": 304.205360}

ENZYMES = {
    # name: (cleave after, not before, or None for N-terminal cleavers)
    "Trypsin": ("KR", "P", "C"),
    "Trypsin/P": ("KR", "", "C"),
    "Lys-C": ("K", "", "C"),
    "Arg-C": ("R", "P", "C"),
    "Chymotrypsin": ("FYWL", "P", "C"),
    "Glu-C": ("E", "", "C"),
    "Asp-N": ("D", "", "N"),
}


def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    acc, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if acc is not None:
                    yield acc, "".join(seq)
                acc, seq = line[1:].split()[0] if len(line) > 1 else "", []
            elif line:
                seq.append(line)
    if acc is not None:
        yield acc, "".join(seq)


def is_decoy(acc):
    a = acc.upper()
    return a.startswith(("DECOY_", "REV_", "XXX_", "DECOY")) or a.endswith(("_REVERSED", "_DECOY", "_REV"))


def digest(seq, enzyme="Trypsin", missed=1, min_len=7, max_len=30):
    after, not_before, side = ENZYMES[enzyme]
    sites = [0]
    for i in range(len(seq) - 1):
        if side == "C" and seq[i] in after and seq[i + 1] not in not_before:
            sites.append(i + 1)
        elif side == "N" and seq[i + 1] in after:
            sites.append(i + 1)
    sites.append(len(seq))
    out = []
    for i in range(len(sites) - 1):
        for j in range(i + 1, min(i + 2 + missed, len(sites))):
            p = seq[sites[i]:sites[j]]
            if min_len <= len(p) <= max_len:
                out.append(p)
    return out


def canon(s, low_res=False):
    s = s.replace("I", "L")
    if low_res:
        s = s.replace("Q", "K")
    return s


@dataclass
class Mods:
    cys: float = 57.021464
    k: float = 0.0          # fixed mass on K (label)
    nterm: float = 0.0      # fixed mass on the peptide N-term (label)

    def residue(self, aa):
        m = AA.get(aa)
        if m is None:
            return None
        if aa == "C":
            m += self.cys
        elif aa == "K":
            m += self.k
        return m


@dataclass
class PeptideIndex:
    name: str
    seqs: list = field(default_factory=list)
    canon: list = field(default_factory=list)
    masses: np.ndarray = None
    decoy: np.ndarray = None   # bool per peptide


def build_index(name, proteins, mods: Mods, enzyme="Trypsin", missed=1, low_res=False,
                max_proteins=None):
    """Digest targets and reversed decoys of a protein list [(acc, seq)]."""
    seen = {}
    for n, (acc, seq) in enumerate(proteins):
        if max_proteins and n >= max_proteins:
            break
        if is_decoy(acc):
            continue
        seq = re.sub(r"[^A-Z]", "", seq.upper())
        for dec, s in ((False, seq), (True, seq[::-1])):
            for p in digest(s, enzyme, missed):
                if p in seen and (seen[p] is False or dec):
                    continue
                seen[p] = dec
    seqs, masses, dec = [], [], []
    for p, d in seen.items():
        m = 0.0
        ok = True
        for aa in p:
            r = mods.residue(aa)
            if r is None:
                ok = False
                break
            m += r
        if not ok:
            continue
        seqs.append(p)
        masses.append(m + H2O + mods.nterm)
        dec.append(d)
    order = np.argsort(masses)
    idx = PeptideIndex(name)
    idx.seqs = [seqs[i] for i in order]
    idx.canon = [canon(seqs[i], low_res) for i in order]
    idx.masses = np.asarray(masses)[order]
    idx.decoy = np.asarray(dec, dtype=bool)[order]
    return idx


def spectrum_tags(scan, frag_tol, mods: Mods, low_res=False, top=60, min_len=3, max_tags=20):
    """Sequence tags: paths through peaks whose spacings are residue masses.
    Returns [(tag, start_mz, end_mz)] where the tag reads from low to high m/z."""
    if scan.mz.size < 5:
        return []
    idx = np.argsort(scan.intensity)[::-1][:top]
    mz = np.sort(scan.mz[idx])
    residues = {}
    for aa in "GASPVTCLNDQKEMHFRYW":
        m = mods.residue(aa)
        residues.setdefault(round(m, 3), set()).add(canon(aa, low_res))
    rmass = np.array(sorted(residues))
    rlett = [sorted(residues[m])[0] for m in rmass]
    edges = defaultdict(list)
    for i in range(mz.size):
        d = mz[i + 1:] - mz[i]
        for j_off in np.nonzero((d > 56) & (d < 187))[0]:
            j = i + 1 + j_off
            k = np.searchsorted(rmass, d[j_off])
            for kk in (k - 1, k):
                if 0 <= kk < rmass.size and abs(rmass[kk] - d[j_off]) <= frag_tol:
                    edges[i].append((j, rlett[kk]))
    # longest paths by DFS from each node (graphs are small)
    best = {}

    def walk(i, tag, start):
        ext = False
        for j, a in edges.get(i, ()):
            if len(tag) < 12:
                ext = True
                walk(j, tag + a, start)
        if not ext and len(tag) >= min_len:
            key = (tag, round(mz[start], 2))
            best[key] = (tag, mz[start], mz[i])

    for i in list(edges):
        walk(i, "", i)
        if len(best) > 400:
            break
    tags = sorted(best.values(), key=lambda t: -len(t[0]))
    return tags[:max_tags]


def _prefix_masses(seq, mods: Mods):
    out = [mods.nterm]
    for aa in seq:
        out.append(out[-1] + (mods.residue(aa) or 0.0))
    return out


def match_spectrum(scan, tags, index: PeptideIndex, prec_ppm, frag_tol, mods: Mods,
                   charges=(2, 3), low_res=False):
    """Return (peptide, is_decoy, precursor_ppm_error) of the first consistent
    candidate, or None."""
    if not tags or scan.precursor_mz <= 0:
        return None
    zs = (scan.charge,) if scan.charge > 0 else charges
    for z in zs:
        mass = (scan.precursor_mz - PROTON) * z
        for iso in (0, 1):
            m = mass - iso * ISO
            tol = m * prec_ppm * 1e-6
            lo = np.searchsorted(index.masses, m - tol)
            hi = np.searchsorted(index.masses, m + tol)
            for c in range(lo, hi):
                cs = index.canon[c]
                for tag, start, _end in tags:
                    for t, ion in ((tag, "b"), (tag[::-1], "y")):
                        pos = cs.find(t)
                        while pos >= 0:
                            pm = _prefix_masses(index.seqs[c], mods)
                            if ion == "b":
                                expect = pm[pos] + PROTON
                            else:
                                expect = (pm[-1] - pm[pos + len(t)]) + H2O + PROTON
                            if abs(expect - start) <= max(frag_tol * 2, 0.02):
                                err = (m - index.masses[c]) / index.masses[c] * 1e6
                                return index.seqs[c], bool(index.decoy[c]), err
                            pos = cs.find(t, pos + 1)
    return None


def search(scans, index, prec_ppm, frag_tol, mods, low_res=False, tag_cache=None):
    """Search spectra against one index. Returns dict with target/decoy counts,
    the matched peptides and precursor errors of target matches."""
    t = d = 0
    hits = {}
    errs = []
    for k, s in enumerate(scans):
        tags = tag_cache[k] if tag_cache is not None else spectrum_tags(s, frag_tol, mods, low_res)
        r = match_spectrum(s, tags, index, prec_ppm, frag_tol, mods, low_res=low_res)
        if r is None:
            continue
        if r[1]:
            d += 1
        else:
            t += 1
            hits[k] = r[0]
            errs.append(r[2])
    return {"target": t, "decoy": d, "net": t - d, "hits": hits, "errors": errs}


def precompute_tags(scans, frag_tol, mods, low_res=False):
    return [spectrum_tags(s, frag_tol, mods, low_res) for s in scans]


def select_scans(run, n=2000):
    ms2 = [s for s in run.ms2 if s.precursor_mz > 0 and s.mz.size >= 10]
    ms2.sort(key=lambda s: -s.tic)
    return ms2[:n]


def load_panel(path):
    """A directory of FASTA files, one per organism (file name = organism)."""
    panel = {}
    for f in sorted(os.listdir(path)):
        if re.search(r"\.(fasta|fa|faa)(\.gz)?$", f, re.I):
            name = re.sub(r"\.(fasta|fa|faa)(\.gz)?$", "", f, flags=re.I)
            panel[name] = os.path.join(path, f)
    return panel


def fdr_net(res):
    """Decoy-corrected matches and a crude FDR."""
    return res["net"], (res["decoy"] / res["target"]) if res["target"] else 1.0


def robust_offset(errs):
    if len(errs) < 20:
        return None
    e = np.asarray(errs)
    med = float(np.median(e))
    mad = float(1.4826 * np.median(np.abs(e - med)))
    return med, mad, len(e)


# --- scoring mini-search -----------------------------------------------------
# Tags alone are too strict on low-resolution or sparse spectra. The mini
# search below scores every candidate in the precursor window by its matched
# b/y ions (a binomial probability score, as in Andromeda/MS Amanda) and
# estimates the FDR from the reversed decoys in the same index.

def _fragments(seq, mods: Mods, max_charge=1):
    pm = np.cumsum([mods.nterm] + [mods.residue(a) or 0.0 for a in seq])
    total = pm[-1] + H2O
    b = pm[1:-1] + PROTON
    y = total - pm[1:-1] + PROTON
    frags = [b, y]
    if max_charge >= 2:
        frags += [(b + PROTON) / 2, (y + PROTON) / 2]
    return np.sort(np.concatenate(frags))


def _binom_score(k, n, p):
    """-log10 P(X >= k) for X ~ Binomial(n, p)."""
    if k <= 0:
        return 0.0
    from math import comb, log10
    tail = 0.0
    for i in range(k, n + 1):
        tail += comb(n, i) * p ** i * (1 - p) ** (n - i)
    return -log10(max(tail, 1e-300))


def prepare_peaks(scan, top_per_100=10):
    """Keep the top peaks per 100 Th window (as most engines do)."""
    if scan.mz.size == 0:
        return scan.mz
    keep = []
    bins = (scan.mz // 100).astype(int)
    for b in np.unique(bins):
        sel = np.nonzero(bins == b)[0]
        keep.extend(sel[np.argsort(scan.intensity[sel])[::-1][:top_per_100]])
    return np.sort(scan.mz[np.array(keep, dtype=int)])


def score_search(scans, index, prec_ppm, frag_tol, frag_ppm, mods: Mods,
                 charges=(2, 3), peaks=None):
    """Best-scoring candidate per spectrum. Returns a list of
    (scan_index, score, is_decoy, peptide, precursor_ppm_error)."""
    out = []
    cache = {}
    for k, s in enumerate(scans):
        pk = peaks[k] if peaks is not None else prepare_peaks(s)
        if pk.size < 5 or s.precursor_mz <= 0:
            continue
        rng = max(1.0, pk[-1] - pk[0])
        best = None
        zs = (s.charge,) if s.charge > 0 else charges
        for z in zs:
            mass = (s.precursor_mz - PROTON) * z
            for iso in (0, 1):
                m = mass - iso * ISO
                tol = m * prec_ppm * 1e-6
                lo = np.searchsorted(index.masses, m - tol)
                hi = np.searchsorted(index.masses, m + tol)
                for c in range(lo, min(hi, lo + 500)):
                    key = (c, z >= 3)
                    th = cache.get(key)
                    if th is None:
                        th = _fragments(index.seqs[c], mods, 2 if z >= 3 else 1)
                        if len(cache) < 200000:
                            cache[key] = th
                    j = np.searchsorted(pk, th)
                    j = np.clip(j, 1, pk.size - 1)
                    near = np.minimum(np.abs(pk[j] - th), np.abs(pk[j - 1] - th))
                    t = th * frag_tol * 1e-6 if frag_ppm else frag_tol
                    hits = int(np.count_nonzero(near <= t))
                    p = min(0.5, pk.size * 2 * (np.mean(t) if frag_ppm else t) / rng)
                    sc = _binom_score(hits, th.size, p)
                    if best is None or sc > best[1]:
                        err = (m - index.masses[c]) / index.masses[c] * 1e6
                        best = (k, sc, bool(index.decoy[c]), index.seqs[c], err)
        if best is not None:
            out.append(best)
    return out


def fdr_filter(psms, fdr=0.01):
    """Target PSMs passing a target-decoy FDR (decoys/targets above the
    threshold). Returns the passing target PSMs and the score threshold."""
    psms = sorted(psms, key=lambda x: -x[1])
    t = d = 0
    best_n, thr = 0, float("inf")
    for i, p in enumerate(psms):
        if p[2]:
            d += 1
        else:
            t += 1
        if t and d / t <= fdr and t > best_n:
            best_n, thr = t, p[1]
    passing = [p for p in psms if p[1] >= thr and not p[2]]
    return passing, thr
