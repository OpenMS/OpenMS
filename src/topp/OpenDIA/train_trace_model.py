#!/usr/bin/env python3
"""Local experimental raw-trace CNN. Outputs are logits, never ID probabilities."""
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Team $
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import torch
from torch import nn

GEOMETRY = [6, 32, 8]
PREPROCESS = "log1p-per-window-max-library-sum-v2"
FIELDS = "sample_id precursor_id group_id is_decoy label rt_apex has_im library_intensities fragment_mz".split()


def protect_inputs(prefix, outputs, extra_inputs=()):
    inputs = [Path(str(prefix) + suffix) for suffix in (".traces.npy", ".samples.tsv", ".geometry.json")]
    inputs.extend(Path(path) for path in extra_inputs)
    outputs = [Path(path) for path in outputs]
    for index, output in enumerate(outputs):
        for other in inputs + outputs[:index]:
            if output.resolve() == other.resolve() or (output.exists() and other.exists() and output.samefile(other)):
                raise ValueError("outputs must be distinct from inputs and each other")


def load_geometry(prefix):
    """Require the exporter's numeric acquisition geometry, not just tensor shape."""
    geometry = json.loads(Path(str(prefix) + ".geometry.json").read_text())
    expected = {"version": 1, "fragments": 6, "rt_bins": 32, "im_bins": 8}
    numeric = ("rt_half_width_seconds", "im_half_width_vssc", "mz_tolerance_da")
    if not isinstance(geometry, dict) or set(geometry) != set(expected) | set(numeric):
        raise ValueError("geometry requires version, dimensions and numeric sampling tolerances")
    if any(type(geometry[k]) is not int or geometry[k] != v for k, v in expected.items()):
        raise ValueError("unsupported geometry version or dimensions")
    if any(type(geometry[k]) not in (int, float) or not np.isfinite(geometry[k]) or geometry[k] <= 0 for k in numeric):
        raise ValueError("geometry sampling tolerances must be finite positive numbers")
    return geometry


class TraceBatches:
    """Read-only disk-backed raw data. Normalize only the requested batch."""
    def __init__(self, raw, geometry):
        self.raw, self.geometry, self.shape = raw, geometry, raw.shape

    def __len__(self):
        return len(self.raw)

    def __getitem__(self, indices):
        x = np.log1p(self.raw[indices])
        x /= np.maximum(x.max(axis=(-3, -2, -1), keepdims=True), 1e-12)
        return torch.from_numpy(x)


def load_dataset(prefix):
    """Validate in bounded tensor batches, without deriving labels from targets."""
    geometry = load_geometry(prefix)
    raw = np.load(str(prefix) + ".traces.npy", allow_pickle=False, mmap_mode="r")
    if raw.dtype.str != "<f4" or raw.ndim != 4 or list(raw.shape[1:]) != GEOMETRY or not raw.flags.c_contiguous:
        raise ValueError("traces must be little-endian C-contiguous float32 [N,6,32,8]")
    for start in range(0, len(raw), 128):
        batch = raw[start:start + 128]
        if not np.isfinite(batch).all() or (batch < 0).any():
            raise ValueError("traces must be finite and nonnegative")
    with open(str(prefix) + ".samples.tsv", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(FIELDS).issubset(reader.fieldnames):
            raise ValueError("samples.tsv missing required columns including group_id")
        rows = list(reader)
    if len(rows) != len(raw) or not rows:
        raise ValueError("metadata and tensor counts must match and be nonempty")
    libraries, masks, labels = [], [], []
    seen = set()
    for row in rows:
        try:
            if not row["sample_id"] or row["sample_id"] in seen or not row["group_id"].strip() or not row["precursor_id"]:
                raise ValueError("unique sample_id and nonempty precursor_id/group_id required")
            seen.add(row["sample_id"])
            label, mask, decoy = int(row["label"]), int(row["has_im"]), int(row["is_decoy"])
            lib = np.asarray([float(v) for v in row["library_intensities"].split(",")], dtype=np.float32)
            mz = np.asarray([float(v) for v in row["fragment_mz"].split(",")])
            if label not in (-1, 0, 1) or mask not in (0, 1) or decoy not in (0, 1):
                raise ValueError("label must be -1/0/1 and masks must be 0/1")
            if lib.shape != (6,) or mz.shape != (6,) or not np.isfinite(lib).all() or not np.isfinite(mz).all():
                raise ValueError("library intensities and fragment m/z require six finite values")
            if decoy and label == 1:
                raise ValueError("decoy positive seeds are forbidden")
            if (lib < 0).any() or (lib[:-1] < lib[1:]).any() or (mz < 0).any() or not np.isfinite(float(row["rt_apex"])):
                raise ValueError("invalid sorted library intensities, m/z or apex")
            padding = mz == 0
            if (lib[padding] != 0).any() or np.any(raw[len(labels), padding] != 0):
                raise ValueError("zero-padded fragment m/z requires zero library and raw intensity")
        except (AttributeError, TypeError, KeyError, ValueError) as error:
            raise ValueError(f"malformed metadata row {len(labels) + 1}: {error}") from error
        libraries.append(lib)
        masks.append(mask)
        labels.append(label)
    # Only raw traces, library intensity and the missing-mobility mask enter the network.
    lib = np.asarray(libraries)
    lib /= np.maximum(lib.max(axis=1, keepdims=True), 1e-12)
    lib /= np.maximum(lib.sum(axis=1, keepdims=True), 1e-12)
    return rows, TraceBatches(raw, geometry), torch.from_numpy(np.column_stack((lib, masks)).astype(np.float32)), np.asarray(labels)


def group_split(groups, seed):
    """Hash groups, not rows, into fixed 60/20/20 partitions."""
    def partition(group):
        bucket = int.from_bytes(hashlib.sha256(f"{seed}\0{group}".encode()).digest()[:8], "big") % 10
        return "train" if bucket < 6 else "validation" if bucket < 8 else "heldout"
    return np.asarray([partition(group) for group in groups])


class TraceCNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv = nn.Sequential(nn.Conv2d(6, 8, 3, padding=1), nn.ReLU(), nn.AvgPool2d(2),
                                  nn.Conv2d(8, 12, 3, padding=1), nn.ReLU(), nn.AvgPool2d(2))
        self.head = nn.Sequential(nn.Linear(12 * 8 * 2 + 7, 24), nn.ReLU(), nn.Linear(24, 1))

    def forward(self, traces, auxiliary):
        return self.head(torch.cat((self.conv(traces).flatten(1), auxiliary), dim=1)).squeeze(1)


def predict(model, x, auxiliary, batch_size=128):
    model.eval()
    with torch.no_grad():
        return torch.cat([model(x[i:i + batch_size], auxiliary[i:i + batch_size])
                          for i in range(0, len(x), batch_size)]).numpy()


def auc(labels, logits):
    """Pairwise rank AUC with half credit for ties, without sklearn."""
    positives, negatives = logits[labels == 1], logits[labels == 0]
    if not len(positives) or not len(negatives):
        return None
    ordered = np.sort(negatives)
    left = np.searchsorted(ordered, positives, side="left")
    right = np.searchsorted(ordered, positives, side="right")
    return float(np.sum((left + right) / 2) / (len(positives) * len(negatives)))


def write_scores(path, rows, logits, splits):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "group_id", "logit", "split"])
        for row, logit, split in zip(rows, logits, splits):
            writer.writerow([row["sample_id"], row["group_id"], float(logit), split])


def train(prefix, checkpoint, output, seed=42, epochs=40, patience=6, batch_size=32):
    protect_inputs(prefix, (checkpoint, str(checkpoint) + ".summary.json", output))
    if not 1 <= epochs <= 100 or not 1 <= patience <= 20 or batch_size < 1:
        raise ValueError("epochs must be 1..100, patience 1..20 and batch_size positive")
    torch.set_num_threads(5)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    rows, x, auxiliary, labels = load_dataset(prefix)
    groups = [row["group_id"] for row in rows]
    splits = group_split(groups, seed)
    indices = {s: np.flatnonzero((splits == s) & (labels >= 0)) for s in ("train", "validation", "heldout")}
    if set(labels[indices["train"]]) != {0, 1}:
        raise ValueError("training partition requires explicit seed positives AND negatives (labels 1 and 0)")
    if not len(indices["validation"]):
        raise ValueError("validation partition requires explicit labels for early stopping")
    model = TraceCNN()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.003)
    loss_fn = nn.BCEWithLogitsLoss()
    y = torch.tensor(labels, dtype=torch.float32)
    best_loss, best_state, stale = float("inf"), None, 0
    for epoch in range(epochs):
        model.train()
        order = indices["train"][torch.randperm(len(indices["train"])).numpy()]
        for start in range(0, len(order), batch_size):
            batch = order[start:start + batch_size]
            optimizer.zero_grad()
            loss = loss_fn(model(x[batch], auxiliary[batch]), y[batch])
            loss.backward()
            optimizer.step()
        model.eval()
        with torch.no_grad():
            val = indices["validation"]
            val_loss = 0.0
            for start in range(0, len(val), batch_size):
                batch = val[start:start + batch_size]
                val_loss += float(loss_fn(model(x[batch], auxiliary[batch]), y[batch])) * len(batch)
            val_loss /= len(val)
        if val_loss < best_loss - 1e-6:
            best_loss, best_state, stale = val_loss, copy.deepcopy(model.state_dict()), 0
        else:
            stale += 1
        if stale >= patience:
            break
    model.load_state_dict(best_state)
    logits = predict(model, x, auxiliary)
    provenance = {s: sorted({g for g, split in zip(groups, splits) if split == s}) for s in indices}
    summary = {"seed": seed, "torch_version": str(torch.__version__), "numpy_version": np.__version__,
               "threads": 5, "epochs_run": epoch + 1, "epoch_cap": epochs, "patience": patience,
               "batch_size": batch_size, "best_validation_loss": best_loss,
               "sample_counts": {s: int(np.sum(splits == s)) for s in indices},
               "labeled_counts": {s: len(v) for s, v in indices.items()},
               "positive_counts": {s: int(np.sum(labels[v] == 1)) for s, v in indices.items()},
               "unknown_count": int(np.sum(labels < 0)),
               "heldout_auc": auc(labels[indices["heldout"]], logits[indices["heldout"]])}
    torch.save({"state_dict": best_state, "geometry": GEOMETRY, "preprocess": PREPROCESS,
                "sampling": x.geometry,
                "source_groups": provenance, "summary": summary}, checkpoint)
    write_scores(output, rows, logits, splits)
    Path(str(checkpoint) + ".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def score(prefix, checkpoint, output, allow_training_overlap=False):
    protect_inputs(prefix, (output,), (checkpoint,))
    torch.set_num_threads(5)
    saved = torch.load(checkpoint, map_location="cpu", weights_only=True)
    if saved["geometry"] != GEOMETRY or saved["preprocess"] != PREPROCESS:
        raise ValueError("checkpoint tensor geometry/preprocessing mismatch")
    if saved.get("sampling") != load_geometry(prefix):
        raise ValueError("checkpoint sampling geometry mismatch")
    rows, x, auxiliary, _ = load_dataset(prefix)
    provenance = {g: s for s, groups in saved["source_groups"].items() for g in groups}
    splits = [provenance.get(row["group_id"], "external") for row in rows]
    if "train" in splits and not allow_training_overlap:
        raise ValueError("scoring dataset overlaps training groups, refused by default")
    model = TraceCNN()
    model.load_state_dict(saved["state_dict"])
    logits = predict(model, x, auxiliary)
    write_scores(output, rows, logits, splits)
    return logits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("train", "score"):
        p = sub.add_parser(command)
        p.add_argument("--prefix", required=True)
        p.add_argument("--checkpoint", required=True)
        p.add_argument("--output", required=True)
        if command == "train":
            p.add_argument("--seed", type=int, default=42)
            p.add_argument("--epochs", type=int, default=40)
            p.add_argument("--patience", type=int, default=6)
            p.add_argument("--batch-size", type=int, default=32)
        else:
            p.add_argument("--allow-training-overlap", action="store_true",
                           help="diagnostic only: overlapping rows remain labeled train, never heldout")
    args = vars(parser.parse_args())
    command = args.pop("command")
    try:
        result = train(**args) if command == "train" else score(**args)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    if command == "train":
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
