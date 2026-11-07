"""Run with python -m unittest discover -s src/topp/OpenDIA/test -p test_prepare_trace_windows.py."""
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Justin Sing $
import csv
import importlib.util
from pathlib import Path
import shutil
import sqlite3
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("prepare", Path(__file__).resolve().parents[1] / "prepare_trace_windows.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class PrepareTraceWindowsTest(unittest.TestCase):
    def test_labels_family_mapping_and_deterministic_sampling(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            library, teacher = root / "lib.pqp", root / "teacher.osw"
            with sqlite3.connect(library) as c:
                c.executescript("CREATE TABLE PRECURSOR(ID,TRAML_ID,DECOY);"
                                "CREATE TABLE PEPTIDE(ID,UNMODIFIED_SEQUENCE);"
                                "CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID,PEPTIDE_ID);")
                for i in range(7):
                    decoy = i >= 3 and i < 6
                    name = f"DECOY_t{i - 3}" if decoy else f"t{i}"
                    c.execute("INSERT INTO PRECURSOR VALUES (?,?,?)", (i, name, int(decoy)))
                    c.execute("INSERT INTO PEPTIDE VALUES (?,?)", (i, "SHUFFLED" if decoy else f"PEPTIDE{i}"))
                    c.execute("INSERT INTO PRECURSOR_PEPTIDE_MAPPING VALUES (?,?)", (i, i))
            shutil.copyfile(library, teacher)
            with sqlite3.connect(teacher) as c:
                c.executescript("CREATE TABLE FEATURE(ID,PRECURSOR_ID,LEFT_WIDTH,RIGHT_WIDTH,EXP_RT);"
                                "CREATE TABLE SCORE_MS2(FEATURE_ID,RANK,QVALUE);")
                for i in range(3):
                    c.execute("INSERT INTO FEATURE VALUES (?,?,9,11,10)", (i, i))
                    c.execute("INSERT INTO SCORE_MS2 VALUES (?,1,0.0001)", (i,))
            candidates = root / "features.tsv"
            with candidates.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["candidate_row", "precursor_id", "is_decoy", "rt_apex", "area"])
                writer.writerows((i, i, int(3 <= i < 6), 10, 100) for i in range(7))
                writer.writerow([7, 0, 0, 99, 10000])  # An off-peak target is unknown, not a positive.
            output, all_windows = root / "sample.tsv", root / "all.tsv"
            result = m.prepare(str(library), str(teacher), str(candidates), str(output), 6, all_windows=str(all_windows))
            self.assertEqual(result["positive_seeds"], 3)
            self.assertEqual(result["decoys"], 3)
            self.assertEqual(result["peptide_families"], 3)
            first = output.read_bytes()
            m.prepare(str(library), str(teacher), str(candidates), str(output), 6)
            self.assertEqual(output.read_bytes(), first)
            with all_windows.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["label"] for row in rows], ["1", "1", "1", "0", "0", "0", "-1", "-1"])
            self.assertEqual(rows[0]["group_id"], rows[3]["group_id"])
            with self.assertRaisesRegex(ValueError, "distinct"):
                m.prepare(str(library), str(teacher), str(candidates), str(candidates), 6)
            with sqlite3.connect(teacher) as c:
                c.execute("UPDATE PRECURSOR SET TRAML_ID='wrong' WHERE ID=0")
            with self.assertRaisesRegex(ValueError, "IDs differ"):
                m.prepare(str(library), str(teacher), str(candidates), str(output), 6)
            with sqlite3.connect(library) as c:
                c.execute("UPDATE PRECURSOR SET TRAML_ID='unpaired' WHERE ID=3")
            with self.assertRaisesRegex(ValueError, "resolve decoy parent"):
                m.prepare(str(library), str(teacher), str(candidates), str(output), 6)


if __name__ == "__main__":
    unittest.main()
