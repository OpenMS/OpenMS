#!/usr/bin/env python3
"""Visualize TIMS-PASEF centroiding before/after for debugging.

Reads three mzML files produced by FileConverter at three centroiding
settings (off / greedy2d / hillbased) and renders side-by-side (m/z, IM)
scatter plots for one selected MS1 frame and optionally one DIA-MS2 spectrum.
Intensity is encoded as marker size and (log) color.

Typical workflow:

    # Convert .d -> mzML at three settings
    BUILD=/path/to/OpenMS-build
    D=/path/to/run.d

    $BUILD/bin/FileConverter -in $D -out /tmp/off.mzML
    $BUILD/bin/FileConverter -in $D -out /tmp/greedy.mzML \\
      -bruker:ms1_centroid_algo greedy2d \\
      -bruker:ms1_centroid_mz_ppm 5 -bruker:ms1_centroid_im_pct 3 \\
      -bruker:ms2_centroid_algo greedy2d \\
      -bruker:dia_ms2_n_neighbors 1
    $BUILD/bin/FileConverter -in $D -out /tmp/hill.mzML \\
      -bruker:ms1_centroid_algo hillbased \\
      -bruker:ms1_centroid_mz_ppm 5 \\
      -bruker:ms2_centroid_algo hillbased \\
      -bruker:ms2_centroid_mz_ppm 20 \\
      -bruker:dia_ms2_n_neighbors 1

    # Render
    python3 tools/scripts/visualize_pasef_centroiding.py \\
        --off /tmp/off.mzML --greedy /tmp/greedy.mzML --hill /tmp/hill.mzML \\
        --out /tmp/centroid_compare.png

Requires pyopenms and matplotlib. The script picks the first non-empty MS1
frame (and first non-empty DIA-MS2 spectrum if --ms2 is passed) by default;
override with --ms1-index / --ms2-index.
"""

import argparse
import sys
from pathlib import Path

# Heavy imports are deferred to render() so that --help works without them.
np = None
plt = None
oms = None


def _import_runtime_deps():
    global np, plt, oms
    try:
        import matplotlib
        matplotlib.use("Agg")  # works headless
        import matplotlib.pyplot as _plt
        import numpy as _np
    except ImportError as e:
        sys.exit(f"missing dependency: {e} (pip install matplotlib numpy)")
    try:
        import pyopenms as _oms
    except ImportError:
        sys.exit("missing pyopenms — build OpenMS-build/pyOpenMS, "
                 "set PYTHONPATH=OpenMS-build/pyOpenMS, or `pip install pyopenms`")
    np, plt, oms = _np, _plt, _oms


def first_nonempty(exp, ms_level, want_im_array=True, index=None):
    """Pick a spectrum at the given MS level. Falls back to first non-empty
    if `index` is None or out of range; returns None if nothing matches."""
    candidates = [i for i in range(exp.size())
                  if exp[i].getMSLevel() == ms_level and exp[i].size() > 0]
    if not candidates:
        return None
    if index is not None and 0 <= index < len(candidates):
        return candidates[index]
    return candidates[0]


_IM_ARRAY_NAMES = {
    "raw inverse reduced ion mobility array",
    "inverse reduced ion mobility array",
}


def extract_peaks(spec):
    """Return (mz, intensity, im) arrays for a spectrum.

    IM is looked up by FloatDataArray name first so the function picks the
    real per-peak IM array even when the spectrum carries extra arrays
    (e.g. the "im lower bound" / "m/z lower bound" arrays attached by
    `-bruker:expose_hill_bounds true`). Falls back to FloatDataArrays[0]
    on older mzMLs without explicit names. Returns None for IM if the
    spectrum has no per-peak IM array (e.g. DDA MS2 with scalar drift_time).
    """
    mz, intensity = spec.get_peaks()
    fda = spec.getFloatDataArrays()
    im = None
    for k in range(len(fda)):
        try:
            name = fda[k].getName()
        except Exception:
            continue
        if name in _IM_ARRAY_NAMES:
            arr = fda[k]
            im = np.fromiter((arr[i] for i in range(arr.size())),
                              count=arr.size(), dtype=float)
            break
    if im is None and len(fda) > 0:
        try:
            arr = fda[0]
            if not arr.getName():
                im = np.fromiter((arr[i] for i in range(arr.size())),
                                  count=arr.size(), dtype=float)
        except Exception:
            im = None
    return np.asarray(mz, dtype=float), np.asarray(intensity, dtype=float), im


def panel(ax, mz, intensity, im, title, mz_window=None):
    """Render one (m/z, IM) scatter panel. Marker size and color encode
    intensity (log-scaled). Returns the scatter artist for colorbar use, or
    None if the panel was empty / fell back to a 1D stem view."""
    if mz.size == 0:
        ax.set_title(f"{title}\n(empty)")
        return None
    if mz_window is not None:
        lo, hi = mz_window
        mask = (mz >= lo) & (mz <= hi)
        mz, intensity = mz[mask], intensity[mask]
        if im is not None:
            im = im[mask]
    if mz.size == 0:
        ax.set_title(f"{title}\n(empty in m/z window)")
        return None

    if im is None:
        # Fallback: plot intensity vs m/z (1D stem-style)
        ax.scatter(mz, intensity, s=4, alpha=0.5, c="#3070a0")
        ax.set_xlabel("m/z")
        ax.set_ylabel("intensity")
        ax.set_title(f"{title}  (no per-peak IM)\n{mz.size} peaks")
        return None

    log_int = np.log10(np.maximum(intensity, 1.0))
    sizes = 2.0 + 30.0 * (log_int - log_int.min()) / max(1e-9, log_int.max() - log_int.min())
    sc = ax.scatter(mz, im, s=sizes, c=log_int, cmap="viridis", alpha=0.6,
                    linewidths=0)
    ax.set_xlabel("m/z")
    ax.set_ylabel("1/K0 (IM)")
    ax.set_title(f"{title}\n{mz.size} peaks")
    return sc


def render(off_path, greedy_path, hill_path, out_path,
           ms1_index=None, ms2_index=None, mz_window=None):
    _import_runtime_deps()
    paths = {"off": off_path, "greedy": greedy_path, "hill": hill_path}
    exps = {}
    for label, p in paths.items():
        if p is None:
            continue
        e = oms.MSExperiment()
        oms.MzMLFile().load(str(p), e)
        exps[label] = e
        print(f"loaded {label}: {e.size()} spectra ({p})", file=sys.stderr)

    if not exps:
        sys.exit("no input mzMLs provided")

    do_ms2 = ms2_index is not None
    n_rows = 2 if do_ms2 else 1
    n_cols = len(exps)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows),
                             squeeze=False)
    fig.suptitle("PASEF centroiding: m/z × IM scatter (intensity → marker size + color)",
                 fontsize=12)

    def _decode_nid(sp):
        n = sp.getNativeID()
        return n.decode() if isinstance(n, bytes) else n

    def _find_by_nid(exp, target_nid, ms_level):
        """Return the index of the spectrum whose nativeID matches target_nid
        at the requested MS level, or None if not found."""
        for i in range(exp.size()):
            if exp[i].getMSLevel() != ms_level:
                continue
            if _decode_nid(exp[i]) == target_nid:
                return i
        return None

    # Pick the anchor frame on the FIRST input so all columns show the same
    # physical (frame, window) — otherwise off/greedy/hill would each show
    # their own "first non-empty" and could be different acquisitions if
    # centroiding dropped a spectrum. Falls back to per-input first_nonempty
    # for inputs that don't contain the anchor.
    anchor_label, anchor_exp = next(iter(exps.items()))
    anchor_ms1_idx = first_nonempty(anchor_exp, 1, index=ms1_index)
    anchor_ms1_nid = (_decode_nid(anchor_exp[anchor_ms1_idx])
                      if anchor_ms1_idx is not None else None)
    if do_ms2:
        anchor_ms2_idx = first_nonempty(anchor_exp, 2, index=ms2_index)
        anchor_ms2_nid = (_decode_nid(anchor_exp[anchor_ms2_idx])
                          if anchor_ms2_idx is not None else None)

    last_sc = None
    for col, (label, exp) in enumerate(exps.items()):
        idx = None
        if anchor_ms1_nid is not None:
            idx = _find_by_nid(exp, anchor_ms1_nid, 1)
        if idx is None:
            idx = first_nonempty(exp, 1, index=ms1_index)
        if idx is None:
            axes[0, col].set_title(f"{label}: no MS1 spectra")
            continue
        sp = exp[idx]
        mz, inten, im = extract_peaks(sp)
        title = f"{label.upper()} MS1 frame={_decode_nid(sp)}"
        sc = panel(axes[0, col], mz, inten, im, title, mz_window=mz_window)
        if sc is not None:
            last_sc = sc

    if do_ms2:
        for col, (label, exp) in enumerate(exps.items()):
            idx = None
            if anchor_ms2_nid is not None:
                idx = _find_by_nid(exp, anchor_ms2_nid, 2)
            if idx is None:
                idx = first_nonempty(exp, 2, index=ms2_index)
            if idx is None:
                axes[1, col].set_title(f"{label}: no MS2 spectra")
                continue
            sp = exp[idx]
            mz, inten, im = extract_peaks(sp)
            native = _decode_nid(sp)
            title = (f"{label.upper()} MS2 ({native[:60]}...)"
                     if len(native) > 60
                     else f"{label.upper()} MS2 {native}")
            sc = panel(axes[1, col], mz, inten, im, title, mz_window=mz_window)
            if sc is not None:
                last_sc = sc

    if last_sc is not None:
        cbar = fig.colorbar(last_sc, ax=axes.ravel().tolist(), pad=0.02,
                            fraction=0.04)
        cbar.set_label("log10(intensity)")

    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    print(f"wrote {out_path}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--off",    type=Path, help="mzML produced with no centroiding")
    ap.add_argument("--greedy", type=Path, help="mzML produced with Greedy2D centroiding")
    ap.add_argument("--hill",   type=Path, help="mzML produced with HillBased centroiding")
    ap.add_argument("--out",    type=Path, default="centroid_compare.png",
                    help="output image path (default: %(default)s)")
    ap.add_argument("--ms1-index", type=int, default=None,
                    help="0-based index among non-empty MS1 spectra (default: first)")
    ap.add_argument("--ms2-index", type=int, default=None,
                    help="0-based index among non-empty MS2 spectra. If set, "
                         "adds a second row of MS2 panels.")
    ap.add_argument("--mz-window", type=str, default=None,
                    help="restrict view to a m/z range, e.g. '600,650'")
    args = ap.parse_args()

    if not (args.off or args.greedy or args.hill):
        ap.error("provide at least one of --off / --greedy / --hill")

    mz_window = None
    if args.mz_window:
        try:
            lo, hi = (float(x) for x in args.mz_window.split(","))
            mz_window = (lo, hi)
        except ValueError:
            ap.error(f"--mz-window must be 'lo,hi', got {args.mz_window!r}")

    render(args.off, args.greedy, args.hill, args.out,
           ms1_index=args.ms1_index, ms2_index=args.ms2_index,
           mz_window=mz_window)


if __name__ == "__main__":
    main()
