#!/usr/bin/env python3
"""Render N MS1 and N MS2 PASEF frames side-by-side across centroiding settings.

Companion to ``visualize_pasef_centroiding.py``: instead of one frame, this
plots many frames in a grid (rows = frames, columns = inputs) so you can
scroll through 10 MS1 and 10 MS2 spectra at a glance to spot regressions
between e.g. off / greedy / hill centroiding.

Typical workflow::

    BUILD=/path/to/OpenMS-build
    D=/path/to/run.d

    $BUILD/bin/FileConverter -in $D -out /tmp/off.mzML
    $BUILD/bin/FileConverter -in $D -out /tmp/hill.mzML \\
        -bruker:ms1_centroid_algo hillbased \\
        -bruker:ms2_centroid_algo hillbased

    python3 tools/scripts/plot_pasef_frames.py \\
        --mzml off=/tmp/off.mzML --mzml hill=/tmp/hill.mzML \\
        --n-ms1 10 --n-ms2 10 \\
        --out-dir /tmp/frames

Outputs ``frames_ms1.png`` and ``frames_ms2.png`` in ``--out-dir``. Frames are
picked at evenly spaced TIC quantiles by default (for ``--n-ms1 10`` the rows
are the 10%, 20%, ..., 100% TIC ranks, last row = highest-TIC frame). Pass
``--pick index`` to fall back to evenly spaced spectrum indices.

Requires pyopenms + matplotlib + numpy. Heavy imports are deferred so
``--help`` works without them.
"""

import argparse
import sys
from pathlib import Path

np = None
plt = None
oms = None

# Bruker parameters whose values matter for hill/centroiding interpretation
# and that should be shown in the figure's params banner. Other (boilerplate)
# bruker params like `calibrate`, `export_mode`, IO flags etc. are omitted.
_BRUKER_PARAMS_TO_SHOW = (
    "ms1_centroid_algo", "ms1_centroid_mz_ppm", "ms1_centroid_im_pct",
    "ms1_centroid_min_hill_length",
    "ms2_centroid_algo", "ms2_centroid_mz_ppm",
    "ms2_centroid_min_hill_length",
    "centroid_valley_factor", "centroid_max_scan_gap",
    "expose_hill_bounds",
    "ms1_n_neighbors", "ms1_min_support",
    "dia_ms2_n_neighbors", "dia_ms2_min_support", "dia_ms2_centroid",
)


def read_bruker_params(path, head_bytes=400_000):
    """Extract `bruker:*` userParams that FileConverter / PeakPickerIM write
    into the mzML's <dataProcessingList>. Reads only the head of the file so
    it's cheap on large mzMLs. Returns dict of short key (without 'bruker:'
    prefix) -> value string. Returns empty dict if no userParams are found."""
    import re
    try:
        with open(path, "rb") as f:
            head = f.read(head_bytes).decode("utf-8", errors="ignore")
    except OSError:
        return {}
    pat = re.compile(
        r'<userParam\s+name="parameter:\s+bruker:([^"]+)"\s+'
        r'type="[^"]*"\s+value="([^"]*)"', re.S)
    return {m.group(1): m.group(2) for m in pat.finditer(head)}


def _import_runtime_deps():
    global np, plt, oms
    try:
        import matplotlib
        matplotlib.use("Agg")
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


def nonempty_indices(exp, ms_level):
    return [i for i in range(exp.size())
            if exp[i].getMSLevel() == ms_level and exp[i].size() > 0]


def pick_evenly(indices, n):
    if not indices or n <= 0:
        return []
    if n >= len(indices):
        return indices
    # n evenly spaced positions across the available indices, inclusive of ends
    step = (len(indices) - 1) / (n - 1) if n > 1 else 0
    return [indices[int(round(k * step))] for k in range(n)]


def pick_by_tic_quantile(exp, indices, n):
    """Pick n frames at evenly spaced TIC quantiles (10%, 20%, ..., 100% for
    n=10). Position k (1-based) maps to the frame ranked at k/n through the
    TIC-ascending order, so the last pick is always the highest-TIC frame."""
    if not indices or n <= 0:
        return []
    tics = []
    for idx in indices:
        _, inten = exp[idx].get_peaks()
        tics.append((float(np.sum(inten)), idx))
    tics.sort(key=lambda x: x[0])  # ascending
    if n >= len(tics):
        return [idx for _, idx in tics]
    m = len(tics)
    picks = []
    for k in range(1, n + 1):
        rank = int(round(k * m / n)) - 1  # 1-based quantile → 0-based index
        rank = max(0, min(m - 1, rank))
        picks.append(tics[rank])
    # Return ordered by quantile (low → high TIC) so the grid reads
    # bottom-percentile → top-percentile top-to-bottom.
    return [idx for _, idx in picks]


def _array_by_name(fda, names):
    """Return the first FloatDataArray whose name matches any of `names`."""
    for k in range(len(fda)):
        try:
            n = fda[k].getName()
        except Exception:
            continue
        if n in names:
            arr = fda[k]
            return np.fromiter((arr[i] for i in range(arr.size())),
                                count=arr.size(), dtype=float)
    return None


def extract_peaks(spec):
    """Return mz, intensity, im (apex), plus optional bounds arrays for hills.

    The IM_apex array is detected by name first ('raw inverse reduced ion
    mobility array' or 'inverse reduced ion mobility array'); falls back to
    FloatDataArrays[0] for older mzMLs without explicit names.
    Hill-bounds arrays ('im lower bound', 'im upper bound', 'm/z lower bound',
    'm/z upper bound') are returned when present, else None.
    """
    mz, intensity = spec.get_peaks()
    fda = spec.getFloatDataArrays()
    im = _array_by_name(fda,
        {"raw inverse reduced ion mobility array",
         "inverse reduced ion mobility array"})
    if im is None and len(fda) > 0:
        arr = fda[0]
        try:
            if not arr.getName():
                im = np.fromiter((arr[i] for i in range(arr.size())),
                                  count=arr.size(), dtype=float)
        except Exception:
            im = None
    bounds = {
        "im_lo": _array_by_name(fda, {"im lower bound"}),
        "im_hi": _array_by_name(fda, {"im upper bound"}),
        "mz_lo": _array_by_name(fda, {"m/z lower bound"}),
        "mz_hi": _array_by_name(fda, {"m/z upper bound"}),
    }
    return (np.asarray(mz, dtype=float),
            np.asarray(intensity, dtype=float),
            im, bounds)


def panel(ax, mz, intensity, im, title,
          mz_window=None, shared_mz=None, shared_im=None, shared_log_int=None,
          mode="datashade", bins=(800, 280), overlay_top=0,
          hill_bounds=None):
    if mz_window is not None:
        lo, hi = mz_window
        mask = (mz >= lo) & (mz <= hi)
        mz, intensity = mz[mask], intensity[mask]
        if im is not None:
            im = im[mask]

    ax.set_title(title, fontsize=8)
    if mz.size == 0:
        ax.text(0.5, 0.5, "empty", transform=ax.transAxes,
                ha="center", va="center", fontsize=7, color="grey")
        if shared_mz is not None:
            ax.set_xlim(shared_mz)
        if shared_im is not None:
            ax.set_ylim(shared_im)
        return None

    if im is None:
        # DDA MS2 path: no per-peak IM, fall back to (m/z, intensity) stems.
        ax.vlines(mz, 0, intensity, linewidth=0.4, color="#3070a0", alpha=0.7)
        ax.set_yscale("log")
        ax.set_xlabel("m/z", fontsize=7)
        ax.set_ylabel("intensity", fontsize=7)
        if shared_mz is not None:
            ax.set_xlim(shared_mz)
        return None

    log_int = np.log10(np.maximum(intensity, 1.0))
    if shared_log_int is not None:
        lo_li, hi_li = shared_log_int
    else:
        lo_li, hi_li = float(log_int.min()), float(log_int.max())

    # Pre-compute top-N picks before any mode-specific reordering so the
    # overlay later references the original arrays consistently.
    top_overlay = None
    if overlay_top > 0:
        k = min(int(overlay_top), intensity.size)
        idx_top = np.argpartition(intensity, -k)[-k:]
        idx_top = idx_top[np.argsort(intensity[idx_top])]
        top_overlay = (mz[idx_top].copy(),
                       im[idx_top].copy(),
                       log_int[idx_top].copy())

    if mode == "datashade":
        # Datashader-style: one rectangular bin per output pixel of the axes,
        # aggregated by max log10(intensity) so a bright peak in a sea of
        # noise still lights up its pixel (sum would average it down).
        # Bin counts come from the figure's per-panel pixel dims so detail
        # scales with --col-width/--row-height/--dpi.
        n_x, n_y = bins
        mz_lo = shared_mz[0] if shared_mz is not None else float(mz.min())
        mz_hi = shared_mz[1] if shared_mz is not None else float(mz.max())
        im_lo_ = shared_im[0] if shared_im is not None else float(im.min())
        im_hi_ = shared_im[1] if shared_im is not None else float(im.max())
        # Guard against degenerate ranges.
        if mz_hi <= mz_lo:
            mz_hi = mz_lo + 1.0
        if im_hi_ <= im_lo_:
            im_hi_ = im_lo_ + 1e-3
        x_idx = np.clip(((mz - mz_lo) / (mz_hi - mz_lo) * n_x).astype(np.int32),
                        0, n_x - 1)
        y_idx = np.clip(((im - im_lo_) / (im_hi_ - im_lo_) * n_y).astype(np.int32),
                        0, n_y - 1)
        grid = np.full((n_y, n_x), np.nan, dtype=np.float32)
        # Write in ascending-intensity order so the brightest peak overwrites
        # weaker neighbors in the same pixel → max-aggregation.
        order = np.argsort(log_int)
        grid[y_idx[order], x_idx[order]] = log_int[order]
        masked = np.ma.masked_invalid(grid)
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad(color="white", alpha=0.0)
        sc = ax.imshow(masked, extent=(mz_lo, mz_hi, im_lo_, im_hi_),
                       origin="lower", aspect="auto",
                       interpolation="nearest", cmap=cmap,
                       vmin=lo_li, vmax=hi_li)
    else:
        # Sort by intensity ascending so high-intensity peaks are drawn LAST
        # (matplotlib scatter z-orders by array position → last = on top).
        order = np.argsort(intensity)
        mz, im = mz[order], im[order]
        log_int = log_int[order]
        rng = max(1e-9, hi_li - lo_li)
        sizes = 1.5 + 22.0 * np.clip((log_int - lo_li) / rng, 0.0, 1.0)
        sc = ax.scatter(mz, im, s=sizes, c=log_int, cmap="viridis",
                        vmin=lo_li, vmax=hi_li, alpha=0.7, linewidths=0)

    # Hill bounding boxes — one per hill of length >= 2 (single-peak hills,
    # which have im_lower == im_upper, are skipped). Independent of the
    # top-N circle overlay above. Requires the mzML to have been produced
    # with `-bruker:expose_hill_bounds true`.
    if (hill_bounds is not None
            and hill_bounds.get("im_lo") is not None
            and hill_bounds.get("im_hi") is not None
            and hill_bounds.get("mz_lo") is not None
            and hill_bounds.get("mz_hi") is not None):
        from matplotlib.patches import Rectangle
        from matplotlib.collections import PatchCollection
        im_lo_arr = hill_bounds["im_lo"]
        im_hi_arr = hill_bounds["im_hi"]
        mz_lo_arr = hill_bounds["mz_lo"]
        mz_hi_arr = hill_bounds["mz_hi"]
        # Multi-peak mask: a hill must span >= 2 IM scans to be drawn.
        # (Length-1 hills have im_lower == im_upper exactly.)
        multi = im_hi_arr > im_lo_arr
        # Padding for hills with measurable IM span but degenerate m/z
        # (a hill that stayed at exactly one m/z across multiple scans),
        # so the box has a finite width on the m/z axis.
        min_w_mz = (shared_mz[1] - shared_mz[0]) * 0.0004 if shared_mz else 0.005
        idxs = np.nonzero(multi)[0]
        if idxs.size:
            rects = []
            for k in idxs:
                x_lo = float(mz_lo_arr[k]); x_hi = float(mz_hi_arr[k])
                y_lo = float(im_lo_arr[k]); y_hi = float(im_hi_arr[k])
                w = x_hi - x_lo
                if w < min_w_mz:
                    cx = 0.5 * (x_lo + x_hi)
                    x_lo = cx - min_w_mz / 2.0
                    w = min_w_mz
                rects.append(Rectangle((x_lo, y_lo), w, y_hi - y_lo))
            # Bright magenta with a thin white halo underneath for visibility
            # against viridis dark and bright backgrounds.
            halo = PatchCollection(rects, facecolor="none", edgecolor="white",
                                   linewidth=0.9, alpha=0.45, zorder=4)
            top = PatchCollection(rects, facecolor="none", edgecolor="#ff1493",
                                  linewidth=0.45, alpha=0.85, zorder=4)
            ax.add_collection(halo)
            ax.add_collection(top)

    # Top-N intensity overlay: large circles colored by intensity, drawn
    # last so they pop above the datashade/scatter background. Uses the
    # shared log10(intensity) colormap so the color matches the background.
    if top_overlay is not None:
        top_mz, top_im, top_log = top_overlay
        ax.scatter(top_mz, top_im, s=14.0, c=top_log, cmap="viridis",
                   vmin=lo_li, vmax=hi_li,
                   edgecolors="black", linewidths=0.3,
                   alpha=0.6, zorder=5)

    ax.set_xlabel("m/z", fontsize=7)
    ax.set_ylabel("1/K0", fontsize=7)
    if shared_mz is not None:
        ax.set_xlim(shared_mz)
    if shared_im is not None:
        ax.set_ylim(shared_im)
    ax.tick_params(axis="both", labelsize=6)
    return sc


def _native_id(spec):
    n = spec.getNativeID()
    if isinstance(n, bytes):
        n = n.decode()
    return n


def build_native_id_lookup(exp, ms_level):
    """Map nativeID -> spectrum index for spectra at this MS level."""
    out = {}
    for i in range(exp.size()):
        sp = exp[i]
        if sp.getMSLevel() != ms_level:
            continue
        out[_native_id(sp)] = i
    return out


def _format_bruker_params(params):
    """Render a small subset of bruker:* params as a compact one-line summary,
    omitting defaults so the banner stays readable."""
    if not params:
        return "(no bruker:* params in mzML metadata)"
    # Only show non-default values for the params we care about.
    defaults = {
        "ms1_centroid_algo":        "off",
        "ms1_centroid_mz_ppm":      "0.0",
        "ms1_centroid_im_pct":      "0.0",
        "ms1_centroid_min_hill_length": "1",
        "ms2_centroid_algo":        "off",
        "ms2_centroid_mz_ppm":      "0.0",
        "ms2_centroid_min_hill_length": "1",
        "centroid_valley_factor":   "1.3",
        "centroid_max_scan_gap":    "0",
        "expose_hill_bounds":       "false",
        "ms1_n_neighbors":          "0",
        "ms1_min_support":          "0",
        "dia_ms2_n_neighbors":      "0",
        "dia_ms2_min_support":      "1",
        "dia_ms2_centroid":         "false",
    }
    parts = []
    for k in _BRUKER_PARAMS_TO_SHOW:
        v = params.get(k)
        if v is None or v == defaults.get(k):
            continue
        parts.append(f"{k}={v}")
    return " ".join(parts) if parts else "(all bruker:* params at defaults)"


def render_level(exps, ms_level, n_frames, out_path,
                 mz_window=None, pick_strategy="tic",
                 mode="datashade", dpi=150,
                 col_width=7.0, row_height=3.0,
                 overlay_top=0, show_hill_bounds=False,
                 input_paths=None):
    """Render a grid: rows = n_frames, cols = len(exps).

    Frame selection is anchored on the first input (picks the row-set there),
    then looked up in every other input by nativeID — robust against inputs
    that have different spectrum counts (e.g. hill centroiding can drop
    empty-after-filtering frames). Cells where the nativeID isn't present
    in a given input are rendered as 'n/a'.
    """
    first_label, first_exp = next(iter(exps.items()))
    pool = nonempty_indices(first_exp, ms_level)
    if pick_strategy == "tic":
        picks = pick_by_tic_quantile(first_exp, pool, n_frames)
        denom = len(picks)
        quantile_labels = ([f"{int(round(100 * k / denom))}%"
                            for k in range(1, denom + 1)]
                           if denom else [])
    else:
        picks = pick_evenly(pool, n_frames)
        quantile_labels = [None] * len(picks)
    if not picks:
        print(f"no non-empty MS{ms_level} spectra in {first_label}; skipping",
              file=sys.stderr)
        return False

    # nativeID of each picked spectrum (from the anchor input).
    pick_native_ids = [_native_id(first_exp[i]) for i in picks]

    # Per-input nativeID -> index lookup (built once per input).
    lookups = {label: build_native_id_lookup(exp, ms_level)
               for label, exp in exps.items()}

    n_rows = len(picks)
    n_cols = len(exps)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(col_width * n_cols, row_height * n_rows),
                             squeeze=False)
    # Per-panel pixel dimensions (approx, accounting for axis margins) →
    # used by datashade mode as bin counts so 1 bin ≈ 1 pixel.
    bins = (max(50, int(col_width * dpi * 0.80)),
            max(50, int(row_height * dpi * 0.70)))
    subtitle = ("TIC quantiles 1/n…n/n (top row = lowest TIC, bottom = highest)"
                if pick_strategy == "tic" else "evenly spaced spectrum indices")
    plot_opts = (f"mode={mode}  pick={pick_strategy}  "
                 f"overlay_top={overlay_top}  "
                 f"show_hill_bounds={show_hill_bounds}"
                 + (f"  mz_window={mz_window[0]:.1f}..{mz_window[1]:.1f}"
                    if mz_window else ""))
    fig.suptitle(f"PASEF MS{ms_level} frames — m/z × IM (intensity → size + color)\n"
                 f"frame selection: {subtitle}\n"
                 f"plot: {plot_opts}",
                 fontsize=11)

    # Conversion-time params banner: read bruker:* userParams from each input
    # mzML's <dataProcessingList> and render a one-liner per input. Anchored
    # at the bottom of the figure (bbox_inches="tight" trims any excess).
    banner_text = None
    if input_paths:
        banner_lines = []
        for label in exps.keys():
            path = input_paths.get(label, "?")
            params = read_bruker_params(path) if path != "?" else {}
            summary = _format_bruker_params(params)
            banner_lines.append(f"  {label}: {summary}")
        banner_text = ("inputs (bruker conversion params, defaults omitted):\n"
                       + "\n".join(banner_lines))

    # Pre-extract all panels and compute *per-row* shared m/z and 1/K0
    # ranges so each row zooms to its own data band — important for DIA-MS2
    # where each precursor window covers only a narrow IM slice. The
    # log10(intensity) range is shared GLOBALLY so the colormap is
    # consistent across the whole figure (one colorbar to rule them all).
    panel_data = {}  # (row, col_label) -> (mz, inten, im, bounds, idx)
    row_mz_range = [None] * n_frames
    row_im_range = [None] * n_frames
    li_lo, li_hi = np.inf, -np.inf
    for r, native_id in enumerate(pick_native_ids):
        row_mz_lo, row_mz_hi = np.inf, -np.inf
        row_im_lo, row_im_hi = np.inf, -np.inf
        row_has_im = False
        for label, exp in exps.items():
            idx = lookups[label].get(native_id)
            if idx is None:
                continue
            mz, inten, im, bounds = extract_peaks(exp[idx])
            if mz_window is not None:
                lo, hi = mz_window
                m = (mz >= lo) & (mz <= hi)
                mz, inten = mz[m], inten[m]
                if im is not None:
                    im = im[m]
                bounds = {k: (v[m] if v is not None else None)
                          for k, v in bounds.items()}
            panel_data[(r, label)] = (mz, inten, im, bounds, idx)
            if mz.size:
                row_mz_lo = min(row_mz_lo, float(mz.min()))
                row_mz_hi = max(row_mz_hi, float(mz.max()))
                li = np.log10(np.maximum(inten, 1.0))
                li_lo = min(li_lo, float(li.min()))
                li_hi = max(li_hi, float(li.max()))
            if im is not None and im.size:
                row_has_im = True
                row_im_lo = min(row_im_lo, float(im.min()))
                row_im_hi = max(row_im_hi, float(im.max()))
        if row_mz_hi > row_mz_lo:
            row_mz_range[r] = (row_mz_lo, row_mz_hi)
        if row_has_im and row_im_hi > row_im_lo:
            row_im_range[r] = (row_im_lo, row_im_hi)

    shared_log_int = (li_lo, li_hi) if li_hi > li_lo else None

    last_sc = None
    for r, anchor_idx in enumerate(picks):
        native_id = pick_native_ids[r]
        qlabel = quantile_labels[r]
        shared_mz = row_mz_range[r]
        shared_im = row_im_range[r]
        for c, (label, exp) in enumerate(exps.items()):
            ax = axes[r, c]
            data = panel_data.get((r, label))
            if data is None:
                ax.text(0.5, 0.5, f"no MS{ms_level} match\nfor nativeID",
                        transform=ax.transAxes,
                        ha="center", va="center", fontsize=7, color="grey")
                ax.set_title(f"{label}  (missing)", fontsize=8)
                if shared_mz is not None:
                    ax.set_xlim(shared_mz)
                if shared_im is not None:
                    ax.set_ylim(shared_im)
                continue
            mz, inten, im, bounds, idx = data
            tic = float(np.sum(inten)) if mz.size else 0.0
            prefix = f"q{qlabel} " if qlabel else ""
            title = (f"{prefix}{label}  idx={idx}  npks={mz.size}  "
                     f"TIC={tic:.2e}\n{native_id[:48]}")
            sc = panel(ax, mz, inten, im, title,
                       mz_window=mz_window,
                       shared_mz=shared_mz, shared_im=shared_im,
                       shared_log_int=shared_log_int,
                       mode=mode, bins=bins,
                       overlay_top=overlay_top,
                       hill_bounds=bounds if show_hill_bounds else None)
            if sc is not None:
                last_sc = sc

    # Vertical breathing room so per-panel title (2 lines) doesn't overlap
    # the m/z tick labels of the row above.
    fig.subplots_adjust(hspace=0.55, wspace=0.22, top=0.96, bottom=0.04)
    if banner_text:
        # Bottom-left, monospace, small font. bbox_inches="tight" pulls it
        # in tightly without leaving large whitespace gaps.
        fig.text(0.005, 0.002, banner_text, ha="left", va="bottom",
                 fontsize=7, family="monospace")

    if last_sc is not None:
        cbar = fig.colorbar(last_sc, ax=axes.ravel().tolist(),
                            pad=0.015, fraction=0.025)
        cbar.set_label("log10(intensity)", fontsize=8)
        cbar.ax.tick_params(labelsize=7)

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}", file=sys.stderr)
    return True


def parse_mzml_arg(s):
    if "=" not in s:
        raise argparse.ArgumentTypeError(
            f"--mzml expects LABEL=PATH, got {s!r}")
    label, path = s.split("=", 1)
    label = label.strip()
    if not label:
        raise argparse.ArgumentTypeError("--mzml LABEL must be non-empty")
    return label, Path(path)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mzml", action="append", type=parse_mzml_arg,
                    metavar="LABEL=PATH", required=True,
                    help="input mzML with a short label (repeatable). "
                         "Order is preserved as columns left→right. "
                         "Example: --mzml off=/tmp/off.mzML --mzml hill=/tmp/hill.mzML")
    ap.add_argument("--n-ms1", type=int, default=10,
                    help="number of MS1 frames to plot (default: %(default)s; "
                         "0 to skip MS1)")
    ap.add_argument("--n-ms2", type=int, default=10,
                    help="number of MS2 frames to plot (default: %(default)s; "
                         "0 to skip MS2)")
    ap.add_argument("--out-dir", type=Path, default=Path("pasef_frames"),
                    help="output directory (default: %(default)s)")
    ap.add_argument("--mz-window", type=str, default=None,
                    help="restrict view to a m/z range, e.g. '600,650'")
    ap.add_argument("--pick", choices=("tic", "index"), default="tic",
                    help="frame selection strategy. 'tic' (default) picks "
                         "evenly spaced TIC quantiles: for n=10 the rows "
                         "correspond to the 10%%, 20%%, ..., 100%% TIC "
                         "ranks (last row = highest-TIC frame). 'index' "
                         "picks evenly spaced spectrum indices.")
    ap.add_argument("--mode", choices=("datashade", "scatter"),
                    default="datashade",
                    help="panel rendering. 'datashade' (default) bins peaks "
                         "into a 2-D histogram at ~one bin per output pixel, "
                         "aggregated by max log10(intensity) so bright peaks "
                         "remain visible amid dense low-intensity background. "
                         "'scatter' draws individual peaks (intensity-sorted "
                         "so high-intensity points are on top).")
    ap.add_argument("--dpi", type=int, default=150,
                    help="output PNG dpi (default: %(default)s)")
    ap.add_argument("--col-width", type=float, default=7.0,
                    help="figure width in inches per column (default: %(default)s)")
    ap.add_argument("--row-height", type=float, default=3.0,
                    help="figure height in inches per row (default: %(default)s)")
    ap.add_argument("--overlay-top", type=int, default=100,
                    help="overlay the top-N intensity peaks per panel as "
                         "outlined circles colored by intensity. 0 disables. "
                         "(default: %(default)s)")
    ap.add_argument("--show-hill-bounds", action="store_true",
                    help="draw hill bounding boxes for every multi-scan hill "
                         "(im_lower < im_upper) on each panel — independent of "
                         "--overlay-top. Requires the mzML to have been produced "
                         "with FileConverter -bruker:expose_hill_bounds true. "
                         "Off by default.")
    args = ap.parse_args()

    _import_runtime_deps()

    mz_window = None
    if args.mz_window:
        try:
            lo, hi = (float(x) for x in args.mz_window.split(","))
            mz_window = (lo, hi)
        except ValueError:
            ap.error(f"--mz-window must be 'lo,hi', got {args.mz_window!r}")

    exps = {}
    input_paths = {}
    for label, path in args.mzml:
        if label in exps:
            ap.error(f"duplicate label: {label!r}")
        if not path.exists():
            ap.error(f"input not found: {path}")
        e = oms.MSExperiment()
        oms.MzMLFile().load(str(path), e)
        exps[label] = e
        input_paths[label] = str(path)
        print(f"loaded {label}: {e.size()} spectra ({path})", file=sys.stderr)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.n_ms1 > 0:
        render_level(exps, ms_level=1, n_frames=args.n_ms1,
                     out_path=args.out_dir / "frames_ms1.png",
                     mz_window=mz_window, pick_strategy=args.pick,
                     mode=args.mode, dpi=args.dpi,
                     col_width=args.col_width, row_height=args.row_height,
                     overlay_top=args.overlay_top,
                     show_hill_bounds=args.show_hill_bounds,
                     input_paths=input_paths)
    if args.n_ms2 > 0:
        render_level(exps, ms_level=2, n_frames=args.n_ms2,
                     out_path=args.out_dir / "frames_ms2.png",
                     mz_window=mz_window, pick_strategy=args.pick,
                     mode=args.mode, dpi=args.dpi,
                     col_width=args.col_width, row_height=args.row_height,
                     overlay_top=args.overlay_top,
                     show_hill_bounds=args.show_hill_bounds,
                     input_paths=input_paths)


if __name__ == "__main__":
    main()
