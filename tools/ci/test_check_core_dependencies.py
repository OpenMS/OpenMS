#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Timo Sachsenberg $

"""Behavioral tests for the dependency boundary, using a small source-tree fixture."""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import check_core_dependencies as dependencies


class IncludeTests(unittest.TestCase):
    def test_comments_strings_conditionals_and_splicing(self):
        source = '''// #include <OpenMS/Comment.h>
/*
#include <OpenMS/Block.h>
*/
const char* example = R"example(
#include <OpenMS/RawString.h>
)example";
#if 0
# include /* comment */ <include/OpenMS/Inactive.h>
#endif
#include "OpenMS/Quoted.h" // trailing comment
#include \\
<OpenMS/Continued.h>
#include HEADER_MACRO
'''
        self.assertEqual(list(dependencies.includes(source)), [
            (9, "include/OpenMS/Inactive.h", False),
            (11, "OpenMS/Quoted.h", True),
            (12, "OpenMS/Continued.h", False),
            (14, "HEADER_MACRO", None),
        ])

    def test_comment_delimiters_in_strings_do_not_hide_includes(self):
        source = 'const char* url = "https://example.org/*";\n#include <OpenMS/Real.h>\n'
        self.assertEqual(list(dependencies.includes(source)), [(2, "OpenMS/Real.h", False)])

    def test_line_comment_continuation(self):
        source = '// continued comment \\\n#include <OpenMS/Hidden.h>\n#include <OpenMS/Real.h>\n'
        self.assertEqual(list(dependencies.includes(source)), [(3, "OpenMS/Real.h", False)])


class BoundaryTests(unittest.TestCase):
    header = "src/openms/include/OpenMS/KERNEL/Peak.h"
    source = "src/openms/source/KERNEL/Peak.cpp"
    extra_source = "src/openms/source/KERNEL/PeakIO.cpp"
    format_header = "src/openms/include/OpenMS/FORMAT/Reader.h"
    format_source = "src/openms/source/FORMAT/Reader.cpp"
    private_header = "src/openms/include/OpenMS/FORMAT/Private.h"

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        files = {
            "src/openms/CMakeLists.txt": "openms_add_library(HEADER_FILES ${OpenMS_sources_h})",
            "src/openms/includes.cmake": """
                include(include/OpenMS/KERNEL/sources.cmake)
                if(WITH_READER)
                  include(include/OpenMS/FORMAT/sources.cmake)
                endif()
            """,
            "src/openms/include/OpenMS/KERNEL/sources.cmake": """
                set(directory include/OpenMS/KERNEL)
                set(sources_list_h Peak.h)
            """,
            "src/openms/include/OpenMS/FORMAT/sources.cmake": """
                set(directory include/OpenMS/FORMAT)
                set(sources_list_h)
                if(WITH_READER)
                  list(APPEND sources_list_h Reader.h)
                endif()
                set(private_headers_list_h Private.h)
            """,
            self.header: '#include <OpenMS/config.h>\n#include <vector>\n',
            self.source: '#include <OpenMS/KERNEL/Peak.h>\n',
            self.format_header: '',
            self.format_source: '#include <OpenMS/FORMAT/Reader.h>\n#include <OpenMS/FORMAT/Private.h>\n',
            self.private_header: '',
            "src/openswathalgo/CMakeLists.txt": "openms_add_library(HEADER_FILES ${OpenSwathAlgoHeaders})",
            "src/openswathalgo/source/OPENSWATHALGO/OpenSwathAlgoFiles.cmake": """
                set(header_directory include/OpenMS/OPENSWATHALGO)
                set(header_algo_list ALGO/Algorithm.h)
                set(header_dataaccess_list DATAACCESS/Data.h)
            """,
            "src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Algorithm.h": '',
            "src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/Data.h": '',
        }
        for name, content in files.items():
            self.write(name, content)
        self.seeds = dependencies.SEEDS
        dependencies.SEEDS = (self.header,)
        self.addCleanup(setattr, dependencies, "SEEDS", self.seeds)
        self.core = {self.header, self.source}
        self.baseline = {"schema_version": 1, "edges": []}

    def write(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    def stop(self):
        return {(edge["source"], edge["target"]) for edge in self.baseline["edges"]}

    def scan(self):
        return dependencies.scan(self.root, self.core, self.stop())

    def cli(self, *extra):
        self.write("core_files.txt", "# comment\n" + "".join(f"{name}\n" for name in sorted(self.core)))
        self.write("baseline.json", json.dumps(self.baseline))
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            status = dependencies.main([
                "--root", str(self.root), "--core", str(self.root / "core_files.txt"),
                "--baseline", str(self.root / "baseline.json"),
                "--report", str(self.root / "report.json"), *extra,
            ])
        self.core = dependencies.read_core(self.root / "core_files.txt")
        return status

    def report(self):
        return json.loads((self.root / "report.json").read_text())

    def test_new_edge_fails_cli_and_exact_baseline_allows_only_that_edge(self):
        self.assertEqual(self.cli(), 0)
        self.write(self.source, '#include <include/OpenMS/FORMAT/Reader.h>\n')
        self.assertEqual(self.cli(), 1)
        report = self.report()
        self.assertEqual(report["new_violations"][0]["target"], self.format_header)
        # Without an exception the derived core would absorb the reader and its implementation.
        self.assertEqual(report["new_core_files"], [self.private_header, self.format_header, self.format_source])
        self.baseline["edges"] = [{
            "rule": "core-to-noncore", "source": self.source,
            "target": self.format_header, "reason": "Legacy reader convenience method.",
        }]
        self.assertEqual(self.cli(), 0)
        self.assertEqual(self.report()["new_core_files"], [])
        # An exception for one source does not exempt its header or directory.
        self.write(self.header, '#include "../FORMAT/Reader.h"\n')
        self.assertEqual(self.cli(), 1)

    def test_removed_edge_requires_baseline_cleanup(self):
        self.baseline["edges"] = [{
            "rule": "core-to-noncore", "source": self.source,
            "target": self.format_header, "reason": "Already removed dependency.",
        }]
        self.assertEqual(self.cli(), 1)
        report = self.report()
        self.assertEqual(len(report["stale_baseline"]), 1)
        self.assertEqual(report["new_violations"], [])

    def test_update_core_regenerates_the_snapshot_from_seeds_and_baseline(self):
        self.write(self.source, '#include <OpenMS/FORMAT/Reader.h>\n')
        self.assertEqual(self.cli("--update-core"), 0)
        self.assertTrue((self.root / "core_files.txt").read_text().startswith("# Generated by"))
        self.assertEqual(self.core, {self.header, self.source, self.format_header,
                                     self.format_source, self.private_header})
        self.assertEqual(self.cli(), 0)
        # Baselining the crossing edge shrinks the derived core; the snapshot must follow.
        self.baseline["edges"] = [{
            "rule": "core-to-noncore", "source": self.source,
            "target": self.format_header, "reason": "Legacy reader convenience method.",
        }]
        self.assertEqual(self.cli(), 1)
        self.assertEqual(self.report()["stale_core_files"],
                         [self.private_header, self.format_header, self.format_source])
        self.assertEqual(self.cli("--update-core"), 0)
        self.assertEqual(self.core, {self.header, self.source})

    def test_missing_snapshot_is_an_error_unless_updating(self):
        self.write("baseline.json", json.dumps(self.baseline))
        arguments = ["--root", str(self.root), "--core", str(self.root / "missing.txt"),
                     "--baseline", str(self.root / "baseline.json")]
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(dependencies.main(arguments), 2)
            self.assertEqual(dependencies.main(arguments + ["--update-core"]), 0)
        self.assertEqual(dependencies.read_core(self.root / "missing.txt"), self.core)

    def test_install_lists_include_optional_headers_but_exclude_private_headers(self):
        self.write(self.format_header, '#include "Private.h"\n')
        report = self.scan()
        self.assertTrue(report["files"][self.format_header]["installed"])
        self.assertFalse(report["files"][self.private_header]["installed"])
        self.assertEqual(report["violations"][0]["rule"], "installed-to-private")

    def test_generated_unresolved_and_macro_includes_stay_distinct(self):
        self.write(self.header, '#include <OpenMS/config.h>\n#include <OpenMS/Missing.h>\n#include HEADER\n')
        report = self.scan()
        self.assertEqual([edge["status"] for edge in report["header_includes"]],
                         ["generated", "unresolved", "macro"])
        self.assertEqual(len(report["violations"]), 2)
        self.assertEqual(self.cli(), 1)

    def test_only_mirrored_translation_units_join_the_core(self):
        self.write(self.extra_source, '#include <OpenMS/KERNEL/Peak.h>\n#include <OpenMS/FORMAT/Reader.h>\n')
        report = self.scan()
        self.assertEqual(report["seed_header_reachability"][self.header], [self.header])
        self.assertEqual(report["violations"], [])
        self.assertIn(self.extra_source, report["files"])
        self.assertEqual(report["core_files"], sorted(self.core))
        self.assertEqual(self.cli(), 0)

    def test_cycles_and_runtime_data_are_separate(self):
        self.write(self.header, '#include <OpenMS/FORMAT/Reader.h>\n')
        self.write(self.format_header, '#include <OpenMS/KERNEL/Peak.h>\n')
        self.write(self.source, '// File::find("ignored.obo");\nFile::find("/CV/psi-ms.obo");\n')
        report = self.scan()
        self.assertEqual(report["cycles"], [sorted([self.header, self.format_header])])
        self.assertEqual(report["runtime_resources"], [{
            "source": self.source, "line": 2, "resource": "/CV/psi-ms.obo",
            "evidence": "literal File::find argument",
        }])
        self.assertEqual(report, self.scan())

    def test_missing_core_file_is_an_error_but_update_drops_it(self):
        (self.root / self.source).unlink()
        self.assertEqual(self.cli(), 2)
        self.assertEqual(self.cli("--update-core"), 0)
        self.assertEqual(self.core, {self.header})

    def test_invalid_baseline_is_an_error(self):
        edge = {"rule": "core-to-noncore", "source": self.source, "target": self.format_header, "reason": ""}
        self.baseline["edges"] = [edge]
        self.assertEqual(self.cli(), 2)
        edge["reason"] = "Legacy dependency."
        self.baseline["edges"].append(edge)
        self.assertEqual(self.cli(), 2)

    def test_computed_cmake_list_is_rejected(self):
        self.write("src/openms/include/OpenMS/KERNEL/sources.cmake", """
            set(directory include/OpenMS/KERNEL)
            set(sources_list_h ${computed})
        """)
        self.assertEqual(self.cli(), 2)


class GraphTests(unittest.TestCase):
    def test_cycles_handle_self_edges_disconnected_nodes_and_long_chains(self):
        graph = {str(i): {str(i + 1)} for i in range(2000)}
        graph.update({"2000": {"1999"}, "self": {"self"}, "isolated": set()})
        self.assertEqual(dependencies.cycles(graph), [["1999", "2000"], ["self"]])


if __name__ == "__main__":
    unittest.main()
