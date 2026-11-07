"""Focused offline tests. Run with Python unittest, numpy and torch only."""
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Team $
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest
from unittest import mock

import numpy as np
import torch

spec = importlib.util.spec_from_file_location("trace_model", Path(__file__).resolve().parents[1] / "train_trace_model.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


GEOMETRY = {"version": 1, "fragments": 6, "rt_bins": 32, "im_bins": 8,
            "rt_half_width_seconds": 15, "im_half_width_vssc": 0.05, "mz_tolerance_da": 0.025}


def fixture(prefix, n=160):
    # Equal totals and identical library/mask inputs. Only the raw spatial shape differs.
    raw = np.zeros((n, 6, 32, 8), dtype="<f4")
    rows = []
    for i in range(n):
        label = i % 2
        if label:
            raw[i, :, 14:18, 4] = 10
        else:
            raw[i, :, 2:4, 4] = 10
            raw[i, :, 28:30, 4] = 10
        rows.append(dict(zip(m.FIELDS, [str(i), f"p{i}", f"family{i // 2}", str(1 - label), str(label),
                                        "100", "0", "6,5,4,3,2,1", "400,500,600,700,800,900"])))
    save(prefix, raw, rows)
    return raw, rows


def save(prefix, raw, rows):
    Path(str(prefix) + ".geometry.json").write_text(json.dumps(GEOMETRY))
    np.save(str(prefix) + ".traces.npy", raw)
    with open(str(prefix) + ".samples.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=m.FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TraceModelTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.prefix = self.root / "data"

    def test_group_split(self):
        groups = [f"g{i // 2}" for i in range(200)]
        splits = m.group_split(groups, 42)
        np.testing.assert_array_equal(splits, m.group_split(groups, 42))
        np.testing.assert_array_equal(splits[::2], splits[1::2])
        self.assertEqual(set(splits), {"train", "validation", "heldout"})
        reverse = m.group_split(groups[::-1], 42)
        np.testing.assert_array_equal(splits, reverse[::-1])

    def test_invalid_contract(self):
        raw, rows = fixture(self.prefix, 4)
        for malformed in (raw[:, :, :-1, :], raw.astype("float64"), raw.astype(">f4"), np.asfortranarray(raw)):
            with self.subTest(shape=malformed.shape, dtype=malformed.dtype):
                save(self.prefix, malformed, rows)
                with self.assertRaises(ValueError):
                    m.load_dataset(self.prefix)
        for field, value in (("label", "2"), ("label", "true"), ("has_im", "2"), ("group_id", ""),
                             ("library_intensities", "1,2"), ("fragment_mz", "nan,2,3,4,5,6")):
            bad = [dict(row) for row in rows]
            bad[0][field] = value
            save(self.prefix, raw, bad)
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                m.load_dataset(self.prefix)
        save(self.prefix, raw, rows[:-1])
        with self.assertRaises(ValueError):
            m.load_dataset(self.prefix)
        save(self.prefix, raw, rows)
        path = Path(str(self.prefix) + ".samples.tsv")
        path.write_text(path.read_text().replace("group_id", "wrong"))
        with self.assertRaises(ValueError):
            m.load_dataset(self.prefix)
        path.unlink()
        with self.assertRaises(FileNotFoundError):
            m.load_dataset(self.prefix)

    def test_explicit_seeds_required(self):
        raw, rows = fixture(self.prefix)
        for label in ("-1", "0"):
            for row in rows:
                row["label"] = label
            save(self.prefix, raw, rows)
            with self.assertRaisesRegex(ValueError, "explicit seed positives"):
                m.train(self.prefix, self.root / "model.pt", self.root / "scores.tsv", epochs=1)

    def test_output_preservation(self):
        fixture(self.prefix, 4)
        source = Path(str(self.prefix) + ".samples.tsv")
        original = source.read_bytes()
        checkpoint = self.root / "model.pt"
        for output in (source, checkpoint):
            with self.assertRaisesRegex(ValueError, "distinct from inputs"):
                m.train(self.prefix, checkpoint, output, epochs=1)
        with self.assertRaisesRegex(ValueError, "distinct from inputs"):
            m.score(self.prefix, checkpoint, source)
        alias = self.root / "alias.tsv"
        alias.hardlink_to(source)
        with self.assertRaisesRegex(ValueError, "distinct from inputs"):
            m.train(self.prefix, checkpoint, alias, epochs=1)
        self.assertEqual(source.read_bytes(), original)

    def test_learns_shape_and_round_trip(self):
        raw, rows = fixture(self.prefix)
        np.testing.assert_array_equal(raw.sum(axis=(1, 2, 3)), np.full(len(raw), 240))
        checkpoint, output = self.root / "model.pt", self.root / "scores.tsv"
        summary = m.train(self.prefix, checkpoint, output, seed=42, epochs=25)
        self.assertGreater(summary["heldout_auc"], 0.95)
        self.assertEqual(summary["threads"], 5)
        with output.open() as handle:
            scores = list(csv.DictReader(handle, delimiter="\t"))
        expected = np.asarray([float(row["logit"]) for row in scores])
        with self.assertRaisesRegex(ValueError, "overlaps training"):
            m.score(self.prefix, checkpoint, self.root / "refused.tsv")
        actual = m.score(self.prefix, checkpoint, self.root / "roundtrip.tsv", allow_training_overlap=True)
        np.testing.assert_allclose(actual, expected, rtol=0, atol=0)
        # Same seed reproduces both split and trained CPU parameters exactly.
        second = m.train(self.prefix, self.root / "second.pt", self.root / "second.tsv", seed=42, epochs=25)
        self.assertEqual(summary, second)
        self.assertEqual(output.read_text(), (self.root / "second.tsv").read_text())
        saved = torch.load(checkpoint, weights_only=True)
        self.assertEqual(saved["geometry"], [6, 32, 8])
        self.assertEqual(saved["preprocess"], m.PREPROCESS)
        self.assertEqual(saved["sampling"], GEOMETRY)
        for changed in (dict(GEOMETRY, im_half_width_vssc=0.06), dict(GEOMETRY, rt_half_width_seconds=20)):
            Path(str(self.prefix) + ".geometry.json").write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, "sampling geometry mismatch"):
                m.score(self.prefix, checkpoint, self.root / "wrong_geometry.tsv", allow_training_overlap=True)
        Path(str(self.prefix) + ".geometry.json").write_text(json.dumps(GEOMETRY))
        sets = [set(v) for v in saved["source_groups"].values()]
        self.assertFalse(sets[0] & sets[1] or sets[0] & sets[2] or sets[1] & sets[2])
        # Identifier, decoy, label, m/z and apex fields cannot affect model inputs.
        for row in rows:
            row.update(sample_id="new" + row["sample_id"], precursor_id="other", group_id="external" + row["group_id"],
                       is_decoy="0", label="-1", rt_apex="999", fragment_mz="1,2,3,4,5,6")
        save(self.prefix, raw, rows)
        external = m.score(self.prefix, checkpoint, self.root / "external.tsv")
        np.testing.assert_array_equal(actual, external)
        self.assertIn("external", (self.root / "external.tsv").read_text())

    def test_geometry_required_and_numeric(self):
        fixture(self.prefix, 4)
        path = Path(str(self.prefix) + ".geometry.json")
        path.write_text(json.dumps(dict(GEOMETRY, rt_half_width_seconds=20)))
        self.assertEqual(m.load_dataset(self.prefix)[1].geometry["rt_half_width_seconds"], 20)
        for field, value in (("version", 2), ("fragments", 5), ("im_bins", 9),
                             ("rt_half_width_seconds", 0), ("im_half_width_vssc", 0),
                             ("mz_tolerance_da", -1), ("im_half_width_vssc", "0.05"),
                             ("im_half_width_vssc", float("nan")), ("mz_tolerance_da", True)):
            path.write_text(json.dumps(dict(GEOMETRY, **{field: value})))
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                m.load_dataset(self.prefix)
        path.write_text(json.dumps({k: v for k, v in GEOMETRY.items() if k != "im_half_width_vssc"}))
        with self.assertRaises(ValueError):
            m.load_dataset(self.prefix)
        path.write_text("not json")
        with self.assertRaises(ValueError):
            m.load_dataset(self.prefix)
        path.unlink()
        with self.assertRaises(FileNotFoundError):
            m.load_dataset(self.prefix)

    def test_padding_and_decoy_positive(self):
        raw, rows = fixture(self.prefix, 4)
        rows[1]["is_decoy"] = "1"
        save(self.prefix, raw, rows)
        with self.assertRaisesRegex(ValueError, "decoy positive"):
            m.load_dataset(self.prefix)
        rows[1]["is_decoy"] = "0"
        rows[0]["fragment_mz"] = "400,500,600,700,800,0"
        save(self.prefix, raw, rows)
        with self.assertRaisesRegex(ValueError, "zero-padded"):
            m.load_dataset(self.prefix)
        rows[0]["library_intensities"] = "6,5,4,3,2,0"
        save(self.prefix, raw, rows)
        with self.assertRaisesRegex(ValueError, "zero-padded"):
            m.load_dataset(self.prefix)
        raw[0, 5] = 0
        save(self.prefix, raw, rows)
        _, x, auxiliary, _ = m.load_dataset(self.prefix)
        self.assertEqual(float(auxiliary[0, 5]), 0)
        self.assertEqual(float(x[0][5].sum()), 0)

    def test_mmap_and_batch_preprocessing(self):
        raw, _ = fixture(self.prefix, 513)
        _, x, auxiliary, _ = m.load_dataset(self.prefix)
        self.assertIsInstance(x.raw, np.memmap)
        self.assertFalse(x.raw.flags.writeable)
        log1p = np.log1p
        sizes = []

        def bounded_log1p(batch):
            sizes.append(len(batch))
            self.assertLessEqual(len(batch), 128)
            return log1p(batch)

        with mock.patch.object(m.np, "log1p", side_effect=bounded_log1p):
            result = m.predict(m.TraceCNN(), x, auxiliary)
        self.assertEqual(result.shape, (513,))
        self.assertEqual(sizes, [128, 128, 128, 128, 1])
        expected = np.log1p(raw[:2])
        expected /= expected.max(axis=(1, 2, 3), keepdims=True)
        np.testing.assert_allclose(x[:2], expected)
        # Memory mapping plus normalization must not alter the exporter file.
        np.testing.assert_array_equal(x.raw[:2], raw[:2])

    def test_cli_train_score(self):
        fixture(self.prefix)
        script = str(Path(m.__file__))
        checkpoint = self.root / "cli.pt"
        common = ["--prefix", str(self.prefix), "--checkpoint", str(checkpoint)]
        trained = subprocess.run([sys.executable, script, "train", *common, "--output",
                                  str(self.root / "cli_train.tsv"), "--epochs", "1"],
                                 capture_output=True, text=True)
        self.assertEqual(trained.returncode, 0, trained.stderr)
        self.assertIn('"heldout_auc"', trained.stdout)
        scored = subprocess.run([sys.executable, script, "score", *common, "--output",
                                 str(self.root / "cli_score.tsv"), "--allow-training-overlap"],
                                capture_output=True, text=True)
        self.assertEqual(scored.returncode, 0, scored.stderr)
        self.assertEqual((self.root / "cli_train.tsv").read_text(), (self.root / "cli_score.tsv").read_text())

    def test_missing_im_and_unknown(self):
        raw, rows = fixture(self.prefix, 4)
        rows[0]["label"] = "-1"
        save(self.prefix, raw, rows)
        _, x, auxiliary, labels = m.load_dataset(self.prefix)
        self.assertEqual(x.shape, (4, 6, 32, 8))
        np.testing.assert_array_equal(auxiliary[:, -1], 0)
        self.assertEqual(labels[0], -1)
        self.assertIsNone(m.auc(np.array([1, 1]), np.array([1., 2.])))
        self.assertEqual(m.auc(np.array([0, 1]), np.array([1., 1.])), 0.5)


if __name__ == "__main__":
    unittest.main()
