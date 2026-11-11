"""Compute label-free, profile-based fragment interference features from the stored six-ion XIC arrays.

This is a DIA-NN-inspired diagnostic implementation, not a claim of bitwise equivalence to DIA-NN. For each candidate
window, select the fragment trace with the highest mean Pearson correlation to the other nonempty traces. Smooth that
reference with a 3-bin moving average. Estimate each trace's reference-shaped scale from reference-squared weighted
intensity. If its correlation to the reference is below 0.8, replace that scale with the minimum trace/reference ratio
near the reference apex. Cap each trace at 1.5 times the scaled reference and retain both original and corrected data.
The saved feature matrix is label-free and is intended for pre-plan diagnostics before any model training.

Usage: python interference_profiles.py raw01|raw02
"""
import sys
import time
import numpy as np
import pandas as pd

ROOT = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
CORE = slice(12, 20)
CORR_FALLBACK = 0.8
CAP = 1.5
EPS = 1e-12
FEATURES = ['ref_quality', 'raw_corr_med', 'area_retained', 'corrected_support',
            'corrected_lib_cos', 'raw_to_corrected_cos', 'corrected_core_fraction']


def smooth3(x):
    """Apply a centered three-bin moving average along the last dimension."""
    xp = np.pad(x, ((0, 0), (0, 0), (1, 1)), mode='edge')
    return (xp[:, :, :-2] + xp[:, :, 1:-1] + xp[:, :, 2:]) / 3.0


def corr_rows(x, y):
    """Pearson correlation between paired rows, returning zero for constant profiles."""
    xc = x - x.mean(axis=1, keepdims=True)
    yc = y - y.mean(axis=1, keepdims=True)
    den = np.sqrt((xc * xc).sum(axis=1) * (yc * yc).sum(axis=1))
    return np.divide((xc * yc).sum(axis=1), den, out=np.zeros_like(den), where=den > EPS)


def correct_batch(x, lib):
    """Return the reference-profile-corrected traces and eight summary features for a batch.

    @param[in] x XIC profiles, shape [candidate, fragment, time], nonnegative.
    @param[in] lib Library fragment intensities, shape [candidate, fragment].
    @return Corrected XIC array and feature matrix, both label-free. The corrected profile is a proposed model input,
    while the summaries include raw coherence and correction strength. Do not score corrected-profile coherence because
    the cap operation itself imposes reference-like shape.
    """
    x = np.maximum(np.asarray(x, dtype=np.float32), 0.0)
    lib = np.maximum(np.asarray(lib, dtype=np.float32), 0.0)
    n, nf, nt = x.shape
    # Extracted traces are sparse in some bins, so compare profiles after one-bin-scale smoothing.
    xs = smooth3(x)
    areas = xs.sum(axis=2)
    active = areas > 0

    # Correlate centered chromatograms and choose the active trace most consistently shaped like its peers.
    centered = xs - xs.mean(axis=2, keepdims=True)
    norm = np.sqrt((centered * centered).sum(axis=2))
    denom = norm[:, :, None] * norm[:, None, :]
    corr = np.divide(centered @ centered.transpose(0, 2, 1), denom,
                     out=np.zeros((n, nf, nf), dtype=np.float32), where=denom > EPS)
    pair = active[:, :, None] & active[:, None, :]
    eye = np.eye(nf, dtype=bool)[None, :, :]
    pair &= ~eye
    counts = pair.sum(axis=2)
    quality = np.divide((corr * pair).sum(axis=2), counts,
                        out=np.full((n, nf), -1.0, dtype=np.float32), where=counts > 0)
    quality[~active] = -np.inf
    ref_idx = np.argmax(quality, axis=1)
    row = np.arange(n)
    ref = xs[row, ref_idx]
    ref_sm = smooth3(ref[:, None, :])[:, 0, :]
    ref_max = ref_sm.max(axis=1)

    # Fit each trace to the consensus reference; poorly correlated traces use a local apex ratio.
    ref2 = ref_sm * ref_sm
    weighted_den = (ref_sm * ref2).sum(axis=1)
    scale = np.divide((xs * ref2[:, None, :]).sum(axis=2), weighted_den[:, None],
                      out=np.zeros((n, nf), dtype=np.float32), where=weighted_den[:, None] > EPS)
    ref_tiled = np.broadcast_to(ref_sm[:, None, :], x.shape)
    trace_corr = corr_rows(xs.reshape(-1, nt), ref_tiled.reshape(-1, nt)).reshape(n, nf)
    apex_mask = ref_sm >= (0.5 * ref_max[:, None])
    local_ratio = np.divide(xs, ref_tiled, out=np.full_like(xs, np.inf), where=ref_tiled > EPS)
    local_ratio = np.where(apex_mask[:, None, :], local_ratio, np.inf)
    fallback = np.min(local_ratio, axis=2)
    scale = np.where(trace_corr < CORR_FALLBACK, fallback, scale)
    scale = np.where(np.isfinite(scale), scale, 0.0)

    cap = CAP * scale[:, :, None] * ref_tiled
    corrected = np.minimum(xs, cap)
    corrected = np.where(active[:, :, None], corrected, 0.0).astype(np.float32)

    raw_area, corr_area = areas.sum(axis=1), corrected.sum(axis=(1, 2))
    retained = np.divide(corr_area, raw_area, out=np.ones(n, dtype=np.float32), where=raw_area > EPS)
    corrected_core = corrected[:, :, CORE].sum(axis=2)
    support = (corrected_core > 0).sum(axis=1).astype(np.float32)
    a = np.sqrt(corrected_core)
    b = np.sqrt(lib)
    lib_cos = np.divide((a * b).sum(axis=1), np.sqrt((a * a).sum(axis=1) * (b * b).sum(axis=1)),
                        out=np.zeros(n, dtype=np.float32), where=(a * a).sum(axis=1) * (b * b).sum(axis=1) > EPS)
    raw_core = x[:, :, CORE].sum(axis=2)
    ar, ac = np.sqrt(raw_core), np.sqrt(corrected_core)
    raw_corr_cos = np.divide((ar * ac).sum(axis=1), np.sqrt((ar * ar).sum(axis=1) * (ac * ac).sum(axis=1)),
                             out=np.ones(n, dtype=np.float32), where=(ar * ar).sum(axis=1) * (ac * ac).sum(axis=1) > EPS)
    core_fraction = np.divide(corrected[:, :, CORE].sum(axis=(1, 2)), corr_area,
                              out=np.zeros(n, dtype=np.float32), where=corr_area > EPS)
    corr_med = np.divide((corr * pair).sum(axis=(1, 2)), pair.sum(axis=(1, 2)),
                         out=np.zeros(n, dtype=np.float32), where=pair.sum(axis=(1, 2)) > 0)
    features = np.column_stack((quality[row, ref_idx], corr_med, retained, support, lib_cos,
                                raw_corr_cos, core_fraction)).astype(np.float32)
    features[~np.isfinite(features)] = 0.0
    return corrected, features


def main(raw):
    t0 = time.time()
    xic = np.load(f'{ROOT}/{raw}/xic.npy', mmap_mode='r')
    lib = np.load(f'{ROOT}/{raw}/lib.npy', mmap_mode='r')
    output = np.lib.format.open_memmap(f'{ROOT}/{raw}/interference_features.npy', mode='w+',
                                       dtype=np.float32, shape=(len(xic), len(FEATURES)))
    # Save corrected profiles as a separate array for a future voxel-model input, retaining the original xic.npy.
    corrected_out = np.lib.format.open_memmap(f'{ROOT}/{raw}/xic_corrected.npy', mode='w+', dtype=np.float32, shape=xic.shape)
    batch_size = 30000
    for start in range(0, len(xic), batch_size):
        end = min(start + batch_size, len(xic))
        corrected, features = correct_batch(xic[start:end], lib[start:end])
        corrected_out[start:end] = corrected
        output[start:end] = features
        if start % 150000 == 0:
            print(raw, f'{end}/{len(xic)} ({time.time() - t0:.0f}s)', flush=True)
    output.flush(); corrected_out.flush()
    pop = pd.read_csv(f'{ROOT}/{raw}/population.tsv', sep='\t')
    dec = pop.is_decoy.values == 1
    a = np.asarray(output)
    print(raw, 'done', f'{time.time() - t0:.1f}s', 'features', FEATURES,
          'decoy medians', np.median(a[dec], axis=0).round(3).tolist(),
          'target medians', np.median(a[~dec], axis=0).round(3).tolist(), flush=True)


if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('raw01', 'raw02'):
        raise SystemExit('Usage: interference_profiles.py raw01|raw02')
    main(sys.argv[1])
