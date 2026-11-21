# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Karnajeet Gosavi $

from pathlib import Path
import tempfile
import unittest

from import_std_timings import read_objects, summarize


class ImportStdTimingsTest(unittest.TestCase):
    """Check that timing reports compare the intended compilation commands."""

    def read(self, text):
        """Parse a temporary Ninja log without leaving files in the source tree."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / ".ninja_log"
            path.write_text(text, encoding="utf-8")
            return read_objects(path)

    def test_last_object_entry_without_scan_or_bmi_double_counting(self):
        """Count a rebuilt object once and exclude non-object build steps."""
        objects = self.read(
            "# ninja log v5\n"
            "0\t100\t0\tsrc\\OpenMS.dir\\File.cpp.obj\ta\n"
            "0\t100\t0\tstd.ixx.ifc\tb\n"
            "0\t100\t0\tFile.cpp.obj.ddi\tc\n"
            "200\t450\t0\tsrc/OpenMS.dir/File.cpp.obj\td\n"
            "0\t1000\t0\tstd.ixx.obj\te\n"
        )
        self.assertEqual(objects, {"src/OpenMS.dir/File.cpp.obj": 0.25, "std.ixx.obj": 1.0})

    def test_separates_converted_other_and_unmatched_objects(self):
        """Keep newly added module work separate from matching source work."""
        path = "src/openms/CMakeFiles/OpenMS.dir/source/SYSTEM/File.cpp.obj"
        result = summarize({path: 1.0, "other.obj": 3.0},
                           {path: 1.5, "other.obj": 5.0, "std.ixx.obj": 4.0})
        self.assertIn("| Matching objects | 2 | 2 | 4.000 | 6.500 | +2.500 |", result)
        self.assertIn("| Converted sources | 1 | 1 | 1.000 | 1.500 | +0.500 |", result)
        self.assertIn("| Other matching objects | 1 | 1 | 3.000 | 5.000 | +2.000 |", result)
        self.assertIn("| All object commands (including module objects) | 2 | 3 | 4.000 | 10.500 | +6.500 |", result)

    def test_rejects_unusable_logs(self):
        """Reject malformed timings and logs without comparable outputs."""
        for text in ("bad row\n", "100\t50\t0\ta.obj\thash\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.read(text)
        with self.assertRaises(ValueError):
            summarize({"a.obj": 1.0}, {"b.obj": 1.0})

    def test_file_suffix_does_not_include_unconverted_file_classes(self):
        """Avoid counting MzMLFile.cpp as the converted File.cpp source."""
        prefix = "src/openms/CMakeFiles/OpenMS.dir/source/"
        objects = {prefix + "SYSTEM/File.cpp.obj": 1.0,
                   prefix + "FORMAT/MzMLFile.cpp.obj": 10.0}
        result = summarize(objects, objects)
        self.assertIn("| Converted sources | 1 | 1 | 1.000 | 1.000 | +0.000 |", result)
        self.assertIn("| Other matching objects | 1 | 1 | 10.000 | 10.000 | +0.000 |", result)


if __name__ == "__main__":
    unittest.main()
