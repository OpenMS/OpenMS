// SVG plots of the ProSE application: axes, histograms, bars, an annotated spectrum mirrored
// against the theoretical one, and scatter plots. No dependencies.

const SVG = "http://www.w3.org/2000/svg";

export function el(name, attributes = {}, text) {
  const e = document.createElementNS(SVG, name);
  for (const [k, v] of Object.entries(attributes)) e.setAttribute(k, v);
  if (text !== undefined) e.textContent = text;
  return e;
}

/// The largest of @p values (Math.max(...values) exceeds the call stack for large arrays).
export function maxOf(values) {
  let max = -Infinity;
  for (const v of values) if (v > max) max = v;
  return max;
}

export function minOf(values) {
  let min = Infinity;
  for (const v of values) if (v < min) min = v;
  return min;
}

function niceTicks(min, max, count = 5) {
  const span = max - min || 1;
  const step = 10 ** Math.floor(Math.log10(span / count));
  const factor = [1, 2, 5, 10].find((f) => span / (step * f) <= count) ?? 10;
  const s = step * factor;
  const ticks = [];
  for (let t = Math.ceil(min / s) * s; t <= max + s * 1e-9; t += s) ticks.push(Math.abs(t) < s * 1e-9 ? 0 : t);
  return ticks;
}

function format(v) {
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e5 || a < 1e-3)) return v.toExponential(0);
  return String(+v.toFixed(4));
}

/// Axes of a plot in @p svg for x in [xMin, xMax] and y in [yMin, yMax]; returns the mapping
/// of data to pixel coordinates and the plot area.
export function axes(svg, { xMin, xMax, yMin = 0, yMax, xLabel = "", yLabel = "", yTicks = true }) {
  const [, , width, height] = svg.getAttribute("viewBox").split(" ").map(Number);
  const m = { left: 54, right: 12, top: 14, bottom: 34 };
  const x = (v) => m.left + ((v - xMin) / (xMax - xMin || 1)) * (width - m.left - m.right);
  const y = (v) => height - m.bottom - ((v - yMin) / (yMax - yMin || 1)) * (height - m.top - m.bottom);
  svg.replaceChildren();
  svg.append(el("line", { class: "axis", x1: m.left, x2: width - m.right, y1: height - m.bottom, y2: height - m.bottom }));
  svg.append(el("line", { class: "axis", x1: m.left, x2: m.left, y1: m.top, y2: height - m.bottom }));
  for (const t of niceTicks(xMin, xMax, Math.max(2, Math.floor(width / 100)))) {
    svg.append(el("text", { x: x(t), y: height - m.bottom + 15, "text-anchor": "middle" }, format(t)));
  }
  if (yTicks) {
    for (const t of niceTicks(yMin, yMax, Math.max(2, Math.floor(height / 50)))) {
      svg.append(el("text", { x: m.left - 6, y: y(t) + 4, "text-anchor": "end" }, format(t)));
    }
  }
  svg.append(el("text", { x: width - m.right, y: height - 3, "text-anchor": "end" }, xLabel));
  svg.append(el("text", { x: 4, y: m.top - 2, "text-anchor": "start" }, yLabel));
  return { x, y, area: { left: m.left, right: width - m.right, top: m.top, bottom: height - m.bottom } };
}

function message(svg, text) {
  const [, , width, height] = svg.getAttribute("viewBox").split(" ").map(Number);
  svg.replaceChildren(el("text", { x: width / 2, y: height / 2, "text-anchor": "middle" }, text));
}

function legend(svg, entries, area) {
  const g = el("g", { class: "legend" });
  entries.forEach(({ name, cls }, i) => {
    const yy = area.top + 6 + i * 16;
    g.append(el("rect", { class: cls, x: area.right - 110, y: yy - 9, width: 10, height: 10 }));
    g.append(el("text", { x: area.right - 95, y: yy }, name));
  });
  svg.append(g);
}

/// Overlaid histograms of the series [{name, cls, values}] (common bins).
export function histogram(svg, series, { bins = 40, xLabel = "", yLabel = "PSMs", range } = {}) {
  const all = series.flatMap((s) => s.values).filter(Number.isFinite);
  if (all.length === 0) {
    message(svg, "no data");
    return;
  }
  let [lo, hi] = range ?? [minOf(all), maxOf(all)];
  if (lo === hi) { lo -= 0.5; hi += 0.5; }
  const w = (hi - lo) / bins;
  const counts = series.map((s) => {
    const c = new Array(bins).fill(0);
    for (const v of s.values) {
      if (!Number.isFinite(v) || v < lo || v > hi) continue;
      c[Math.min(bins - 1, Math.floor((v - lo) / w))]++;
    }
    return c;
  });
  const { x, y, area } = axes(svg, { xMin: lo, xMax: hi, yMax: maxOf(counts.flat()) * 1.05, xLabel, yLabel });
  series.forEach((s, i) => {
    counts[i].forEach((n, b) => {
      if (n === 0) return;
      svg.append(el("rect", { class: s.cls, x: x(lo + b * w), y: y(n), width: Math.max(1, x(lo + (b + 1) * w) - x(lo + b * w) - 0.5), height: y(0) - y(n) }));
    });
  });
  if (series.length > 1) legend(svg, series.map(({ name, cls }) => ({ name: `${name} (${series.find((s) => s.name === name).values.length})`, cls })), area);
}

/// Bars at positions x with heights; items [{x, height, cls, label}] (labels of the highest are drawn).
export function bars(svg, items, { xLabel = "", yLabel = "", labels = 10 } = {}) {
  if (items.length === 0) {
    message(svg, "no data");
    return;
  }
  const xs = items.map((i) => i.x);
  const pad = (maxOf(xs) - minOf(xs)) * 0.03 || 1;
  const { x, y } = axes(svg, { xMin: minOf(xs) - pad, xMax: maxOf(xs) + pad, yMax: maxOf(items.map((i) => i.height)) * 1.18, xLabel, yLabel });
  for (const i of items) {
    svg.append(el("rect", { class: `bar ${i.cls ?? ""}`, x: x(i.x) - 1.5, y: y(i.height), width: 3, height: y(0) - y(i.height) }));
  }
  const placed = [];
  for (const i of [...items].sort((a, b) => b.height - a.height).slice(0, labels)) {
    if (!i.label) continue;
    const cx = x(i.x), half = i.label.length * 3;
    if (placed.some(([l, r]) => cx + half > l && cx - half < r)) continue;
    placed.push([cx - half, cx + half]);
    svg.append(el("text", { x: cx, y: y(i.height) - 4, "text-anchor": "middle" }, i.label));
  }
}

/// Scatter plot of points [{x, y, cls, title}], with a dashed line at y = 0.
export function scatter(svg, points, { xLabel = "", yLabel = "", symmetric = true } = {}) {
  if (points.length === 0) {
    message(svg, "no matched ions");
    return;
  }
  const xs = points.map((p) => p.x), ys = points.map((p) => p.y);
  const pad = (maxOf(xs) - minOf(xs)) * 0.05 || 1;
  const yAbs = Math.max(maxOf(ys.map(Math.abs)) * 1.2, 1e-6);
  const { x, y } = axes(svg, {
    xMin: minOf(xs) - pad, xMax: maxOf(xs) + pad,
    yMin: symmetric ? -yAbs : minOf(ys), yMax: symmetric ? yAbs : maxOf(ys), xLabel, yLabel,
  });
  const [xMin, xMax] = [minOf(xs) - pad, maxOf(xs) + pad];
  svg.append(el("line", { class: "zero", x1: x(xMin), x2: x(xMax), y1: y(0), y2: y(0) }));
  for (const p of points) {
    const c = el("circle", { class: `point ${p.cls ?? ""}`, cx: x(p.x), cy: y(p.y), r: 3.5 });
    if (p.title) c.append(el("title", {}, p.title));
    svg.append(c);
  }
}

/// The ion type (a, b, c, x, y, z) of an ion name such as "y7++" or "b3-H2O1+".
export const ionType = (name) => (/^[abcxyz]/.test(name) ? name[0] : "");

/// The measured spectrum (relative intensity, up) with the matched ions in color and labelled,
/// and the theoretical spectrum mirrored below.
export function mirrorSpectrum(svg, peaks, annotation) {
  const { mz, intensity } = peaks;
  if (mz.length === 0) {
    message(svg, "empty spectrum");
    return;
  }
  const theo = annotation?.theoretical ?? { mz: [], ion: [] };
  const allMZ = [...mz, ...theo.mz];
  const pad = (maxOf(allMZ) - minOf(allMZ)) * 0.03 || 1;
  const xMin = minOf(allMZ) - pad, xMax = maxOf(allMZ) + pad;
  const top = maxOf(intensity) || 1;
  const { x, y, area } = axes(svg, { xMin, xMax, yMin: -100, yMax: 118, xLabel: "m/z", yLabel: "relative intensity (%)" });
  svg.append(el("line", { class: "axis", x1: area.left, x2: area.right, y1: y(0), y2: y(0) }));

  for (let i = 0; i < mz.length; ++i) {
    svg.append(el("line", { class: "peak", x1: x(mz[i]), x2: x(mz[i]), y1: y(0), y2: y((intensity[i] / top) * 100) }));
  }
  const matches = annotation?.matches ?? [];
  for (const m of matches) {
    const line = el("line", { class: `peak ${ionType(m.ion)}`, x1: x(m.mz), x2: x(m.mz), y1: y(0), y2: y((m.intensity / top) * 100) });
    line.append(el("title", {}, `${m.ion} m/z ${m.mz.toFixed(4)} (${m.error.toFixed(2)} ${annotation.unit})`));
    svg.append(line);
  }
  // labels: the most intense matches first, left out where they would overlap
  const placed = [];
  for (const m of [...matches].sort((a, b) => b.intensity - a.intensity)) {
    const cx = x(m.mz), cy = y((m.intensity / top) * 100) - 4, half = m.ion.length * 3.1;
    if (placed.some(([l, r, t]) => cx + half > l && cx - half < r && Math.abs(cy - t) < 12)) continue;
    placed.push([cx - half, cx + half, cy]);
    svg.append(el("text", { class: `ion-label ${ionType(m.ion)}`, x: cx, y: cy, "text-anchor": "middle" }, m.ion));
  }
  // theoretical ions, mirrored (all at -70 %, matched ones in color)
  const matchedTheo = new Set(matches.map((m) => m.theoreticalMZ));
  for (let i = 0; i < theo.mz.length; ++i) {
    const cls = matchedTheo.has(theo.mz[i]) ? `peak theoretical ${ionType(theo.ion[i])}` : "peak theoretical";
    const line = el("line", { class: cls, x1: x(theo.mz[i]), x2: x(theo.mz[i]), y1: y(0), y2: y(-70) });
    line.append(el("title", {}, `${theo.ion[i]} m/z ${theo.mz[i].toFixed(4)}`));
    svg.append(line);
  }
  svg.append(el("text", { x: area.left + 6, y: y(-92) }, "theoretical"));
}
