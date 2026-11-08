"""Augmented peak-boundary training examples from clean calibration anchors (step B data; numpy only).

Labels come only from OpenDIACalibrate -out_anchor_peaks (simple consensus picker on raw-scan traces of target anchors).
No decoys, search candidates or external tools are used. Each example is a 6 x 32 window on the same geometry as the
CNN tensor (+-15 s around a centre, 32 bins), with per-bin targets: inside-peak mask and start/end heatmaps.

Augmentations (each sampled per example):
  - RT jitter: the window centre is the picked apex shifted by up to +-6 s (the model must not assume a centred apex);
  - peak scale 10^U(-3, 0) relative to the anchor, with Poisson-like noise sqrt(I); the background below stays at its real
    level, so the scale sets the signal-to-background ratio (recorded per example as log10 SNR);
  - fragment dropout: 0-2 fragments zeroed (missing transitions);
  - background: the peak is pasted onto an off-peak segment of the same anchor's own raw trace (real interference at
    those fragment m/z), taken more than 30 s from the apex, at its native intensity;
  - interference: with p = 0.3 a single fragment gets a foreign Gaussian at another RT;
  - neighbour: with p = 0.2 a second, co-eluting-shaped peak of all fragments 1.5-4 widths away (valley case);
  - negatives: with p = 0.2 the window is background only (no peak; all-zero targets).
The fold of each anchor is the CNN outer fold so a model trained on folds != k never saw fold k's anchors."""
import argparse, hashlib
import numpy as np, pandas as pd

BINS, HALF = 32, 15.0


def fold_of(precursor_id, folds, seed):
    # Deterministic hash split by precursor (the CNN uses its own hash on peptide families; this is a stand-in until wired natively).
    h = hashlib.sha256(f'outer\t{precursor_id}\t{seed}'.encode()).digest()
    return int.from_bytes(h[:8], 'little') % folds


def bin_trace(rt, inten, centre):
    """Sum raw-scan points into the CNN's 32 bins over centre +- 15 s (same rule as TraceCollector::rtBin_)."""
    out = np.zeros((inten.shape[1], BINS), np.float32)
    keep = np.abs(rt - centre) <= HALF
    b = np.minimum(BINS - 1, np.maximum(0.0, (rt[keep] - centre + HALF) / (2 * HALF)) * BINS).astype(int)
    np.add.at(out.T, b, inten[keep])
    return out


def bin_mask(rt, centre):
    """Bins that received at least one raw scan (TraceCollector's ms2_mask)."""
    keep = np.abs(rt - centre) <= HALF
    b = np.minimum(BINS - 1, np.maximum(0.0, (rt[keep] - centre + HALF) / (2 * HALF)) * BINS).astype(int)
    m = np.zeros(BINS, np.int8)
    m[b] = 1
    return m


def bin_of(t, centre):
    return int(np.clip((t - centre + HALF) / (2 * HALF) * BINS, 0, BINS - 1))


def targets(left, right, centre):
    inside = np.zeros(BINS, np.float32)
    start = np.zeros(BINS, np.float32)
    end = np.zeros(BINS, np.float32)
    lo, hi = bin_of(left, centre), bin_of(right, centre)
    inside[lo:hi + 1] = 1
    start[lo], end[hi] = 1, 1
    return inside, start, end


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--peaks', required=True)
    ap.add_argument('--traces', required=True)
    ap.add_argument('--out', required=True, help='.npz')
    ap.add_argument('--per_anchor', type=int, default=64)
    ap.add_argument('--folds', type=int, default=3)
    ap.add_argument('--seed', type=int, default=42)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    peaks = pd.read_csv(a.peaks, sep='\t')
    peaks = peaks[peaks.clean == 1]
    tr = pd.read_csv(a.traces, sep='\t')
    groups = {k: g.sort_values('scan') for k, g in tr.groupby('precursor_id')}
    cols = [f'i{f}' for f in range(6)]
    X, Y, M, meta = [], [], [], []
    for p in peaks.itertuples():
        g = groups[p.precursor_id]
        rt = g.rt.to_numpy()
        inten = g[cols].to_numpy(np.float64)
        # The whole anchor slice is the signal (real tails and neighbours), not a copy truncated at the labelled
        # boundaries: truncation would leave a step at the label edge that a model could learn instead of the peak.
        peak = inten
        width = p.right_rt - p.left_rt
        far = np.abs(rt - p.apex_rt) >= max(3 * width, 2 * HALF)
        bg_rt = rt[far]
        fold = fold_of(p.precursor_id, a.folds, a.seed)
        apex_height = peak.sum(1).max()
        for _ in range(a.per_anchor):
            negative = rng.random() < 0.2
            shift = rng.uniform(-6, 6)
            centre = p.apex_rt + shift
            scale = 10 ** rng.uniform(-3, 0)
            sig = bin_trace(rt, peak, centre) * scale if not negative else np.zeros((6, BINS), np.float32)
            if rng.random() < 0.2 and not negative:
                # Neighbour: the same peak shape displaced by 1.5-4 widths at 0.3-1x height (a second, co-eluting species).
                d = rng.choice([-1, 1]) * rng.uniform(1.5, 4) * max(width, 2.0)
                sig += bin_trace(rt + d, peak, centre) * scale * rng.uniform(0.3, 1.0) * rng.permutation(np.linspace(0.3, 1.0, 6))[:, None].astype(np.float32)
            bg = np.zeros((6, BINS), np.float32)
            if len(bg_rt) > 20:
                c2 = rng.choice(bg_rt)
                if abs(c2 - p.apex_rt) > 2 * HALF:
                    bg = bin_trace(rt, inten, c2)
            x = sig + bg
            # SNR: peak apex over the background consensus median (+1 count) in this window.
            apex_bin = bin_of(p.apex_rt, centre)
            snr = (sig.sum(0)[max(0, apex_bin - 1):apex_bin + 2].max() + 1.0) / (np.median(bg.sum(0)) + 1.0)
            if rng.random() < 0.3:
                f = rng.integers(6)
                mu = rng.uniform(-HALF, HALF)
                t = np.linspace(-HALF, HALF, BINS)
                x[f] += (apex_height * scale * rng.uniform(0.2, 1.5) * np.exp(-0.5 * ((t - mu) / rng.uniform(1, 3)) ** 2)).astype(np.float32)
            drop = rng.integers(0, 3)
            if drop:
                x[rng.choice(6, drop, replace=False)] = 0
            x = np.maximum(0, x + rng.normal(0, 1, x.shape) * np.sqrt(np.maximum(x, 0) + 1e-6)).astype(np.float32)
            if negative:
                y = np.zeros((3, BINS), np.float32)
            else:
                y = np.stack(targets(p.left_rt, p.right_rt, centre))
            X.append(x)
            Y.append(y)
            M.append(bin_mask(rt, centre))
            meta.append((p.precursor_id, fold, int(negative), shift, np.log10(scale), np.log10(snr)))
    m = np.array(meta)
    np.savez_compressed(a.out, x=np.stack(X), y=np.stack(Y), scan_mask=np.stack(M), precursor_id=m[:, 0].astype(np.int64), fold=m[:, 1].astype(np.int8),
                        negative=m[:, 2].astype(np.int8), shift=m[:, 3].astype(np.float32), log_scale=m[:, 4].astype(np.float32),
                        log_snr=m[:, 5].astype(np.float32))
    print(f'{len(X)} examples from {len(peaks)} clean anchors; negatives {int(m[:, 2].sum())}; folds {np.bincount(m[:, 1].astype(int)).tolist()}')


if __name__ == '__main__':
    main()
