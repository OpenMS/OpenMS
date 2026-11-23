"""Synthetic fixtures only. Run: python3 -m unittest discover -s THIS_DIRECTORY -v"""
import json
from pathlib import Path
import tempfile
import unittest

from eval_quant import describe, evaluate, read_run


class QuantTests(unittest.TestCase):
    def test_scaling_missingness_and_zero(self):
        result = describe([{"a": 10., "b": 100., "z": None},
                           {"a": 20., "b": 200.}, {"a": 40., "b": 400.}], ["1", "2", "3"])
        a, b, z = result["features"]
        self.assertAlmostEqual(a["cv_raw"], 0.6546536707079771)
        self.assertAlmostEqual(a["cv_median_normalized"], 0.)
        self.assertEqual(result["positive_complete"], 2)
        self.assertEqual(result["positive_completeness"]["0"], 1)
        self.assertIsNone(z["cv_raw"])
        self.assertEqual(z["raw"], [None, None, None])
        self.assertAlmostEqual(result["pairs"][0]["log2_pearson"], 1.)
        json.dumps(result, allow_nan=False)

    def test_empty_constant_and_single_seed(self):
        self.assertIsNone(describe([{}, {}, {}], ["a", "b", "c"])["pairs"][0]["jaccard"])
        self.assertIsNone(describe([{"x": 1}, {"x": 2}], ["a", "b"])["pairs"][0]["log2_pearson"])
        self.assertIsNone(describe([{"x": 1}], ["s0"])["features"][0]["cv_raw"])
        no_common = describe([{"x": 1}, {"y": 2}], ["a", "b"])
        self.assertEqual(no_common["normalization_log2_offsets"], [None, None])

    def test_population_switches_and_cutoff_grid(self):
        from test_prepare_library import fixture
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture(root / "library.pqp")
            (root / "run.tsv").write_text("run\tprecursor_id\tis_decoy\tsequence\tcharge\tarea\tq_value\n"
                "r\t1\t0\tPEP\t2\t10\t0.004\n"
                "r\t4\t0\tENT\t2\t100\t0.009\n"
                "r\t3\t0\tPEP\t3\t20\t0.015\n")
            manifest = root / "manifest.tsv"
            manifest.write_text("raw\tseed\tpath\n" + "".join(
                f"{r}\t0\trun.tsv\n" for r in ("01", "02", "03")))
            for cutoff, all_n, human_n in ((.005, 1, 1), (.01, 2, 1), (.02, 3, 2)):
                all_result = evaluate(manifest, cutoff, root / "library.pqp", "all")
                human_result = evaluate(manifest, cutoff, root / "library.pqp", "human")
                self.assertEqual(all_result["across_raws"]["precursor"]["0"]["union"], all_n)
                self.assertEqual(human_result["across_raws"]["precursor"]["0"]["union"], human_n)
                self.assertEqual(human_result["population"], "human")
                self.assertEqual(human_result["audit"][0]["population_excluded"], all_n-human_n)
            with self.assertRaises(ValueError):
                evaluate(manifest, .01, population="human")
            with open(root / "run.tsv", "a") as handle:
                handle.write("r\t999\t0\tUNKNOWN\t2\t1\t0\n")
            with self.assertRaisesRegex(ValueError, "not a target"):
                evaluate(manifest, .01, root / "library.pqp", "all")

    def test_schema_filters_manifest_and_aggregation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            header = "run\tprecursor_id\tis_decoy\tsequence\tcharge\tarea\tq_value\n"
            rows = ("r\ta\t0\tPEP\t2\t10\t0.01\n"
                    "r\tb\t0\tPEP\t3\t20\t0.001\n"
                    "r\tz\t0\tZERO\t2\t0\t0.001\n"
                    "r\tn\t0\tNAN\t2\tnan\t0.001\n"
                    "r\tf\t0\tFAIL\t2\t50\t0.02\n"
                    "r\td\t1\tDECOY\t2\t100\t0\n")
            (root / "run.tsv").write_text(header + rows)
            values, counts = read_run(root / "run.tsv", .01)
            self.assertEqual(set(values), {"a", "b", "z", "n"})
            self.assertEqual(counts["nonpositive"], 1)
            self.assertEqual(counts["nonfinite"], 1)
            manifest = root / "manifest.tsv"
            manifest.write_text("raw\tseed\tpath\n" + "".join(
                f"{r}\t{s}\trun.tsv\n" for r in ("01", "02", "03") for s in ("0", "1")))
            result = evaluate(manifest, .01)
            peptide = result["across_raws"]["peptide"]["0"]["features"]
            self.assertEqual(next(x for x in peptide if x["feature"] == "PEP")["raw"], [30., 30., 30.])
            self.assertEqual(len(result["within_raw_seeds"]["precursor"]["01"]["labels"]), 2)
            json.dumps(result, allow_nan=False)
            manifest.write_text("raw\tseed\tpath\n01\t0\trun.tsv\n")
            with self.assertRaises(ValueError):
                evaluate(manifest, .01)
            (root / "run.tsv").write_text(header + rows + "r\ta\t0\tPEP\t2\t1\t0\n")
            with self.assertRaises(ValueError):
                read_run(root / "run.tsv", .01)
            (root / "run.tsv").write_text("precursor_id\tarea\na\t1\n")
            with self.assertRaises(ValueError):
                read_run(root / "run.tsv", .01)


if __name__ == "__main__":
    unittest.main()
