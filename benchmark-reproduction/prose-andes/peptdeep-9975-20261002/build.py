"""Build and class-test each source tree with the native ProSE harness.

Algorithms (ProSEAlgorithm, FragmentIndex, HyperScore, TheoreticalSpectrumGenerator and,
when present, FragmentIonLikelihoodModel) are compiled from sources/<build>; all other
OpenMS code comes from the pinned pyOpenMS wheel's libOpenMS (see README). Test fixtures
are read from the same source tree that is compiled.
"""

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT.parent / "prose-andes-reproduction"
ENV = dict(
    os.environ,
    OMP_NUM_THREADS="2",
    OPENMS_DATA_PATH=str(BUNDLE / "nightly/pyopenms/share/OpenMS"),
)

TEST_CONFIG = """#ifndef OPENMS_TEST_CONFIG_H
#define OPENMS_TEST_CONFIG_H
#define OPENMS_GET_TEST_DATA_PATH(filename) (std::string("{data}/") + filename).c_str()
#define OPENMS_GET_TEST_DATA_PATH_MESSAGE(prefix,filename,suffix) (prefix + std::string("{data}/") + filename + suffix).c_str()
#endif // OPENMS_TEST_CONFIG_H
"""


def git(source, *args):
    return subprocess.check_output(["git", *args], cwd=source, text=True).strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("builds", nargs="+")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--skip-tests", action="store_true")
    ap.add_argument("--harness", default="harness", help="harness directory (harness_onnx: with ONNX Runtime and PeptDeep)")
    args = ap.parse_args()
    record_path = ROOT / "build_records.json"
    for build in args.builds:
        source = ROOT / "sources" / build
        dest = ROOT / "builds" / build
        cfg = dest / "testcfg" / "OpenMS"
        cfg.mkdir(parents=True, exist_ok=True)
        (cfg / "test_config.h").write_text(
            TEST_CONFIG.replace("{data}", str(source / "src/tests/class_tests/openms/data"))
        )
        steps = [
            ["cmake", "-S", ROOT / args.harness, "-B", dest, "-DCMAKE_BUILD_TYPE=Release",
             "-DCMAKE_CXX_COMPILER=g++-13", "-DOPENMS=" + str(source), "-DTESTCFG=" + str(dest / "testcfg")],
            ["cmake", "--build", dest, f"-j{args.jobs}"],
        ]
        if not args.skip_tests:
            steps.append(["ctest", "--test-dir", dest, "--output-on-failure", "-j2"])
        with (dest / "build.log").open("w") as log:
            for cmd in steps:
                print(build, " ".join(map(str, cmd)), flush=True)
                subprocess.run(list(map(str, cmd)), env=ENV, stdout=log, stderr=log, check=True)
        assert not git(source, "status", "--porcelain", "--untracked-files=no"), "dirty source tree: " + build
        records = json.loads(record_path.read_text()) if record_path.exists() else {}
        records[build] = {
            "source_commit": git(source, "rev-parse", "HEAD"),
            "source_tree": git(source, "rev-parse", "HEAD^{tree}"),
            "binary_sha256": hashlib.sha256((dest / "prose").read_bytes()).hexdigest(),
            "tests_passed": not args.skip_tests,
            "harness": args.harness,
        }
        record_path.write_text(json.dumps(records, indent=2))
        print("Build ok:", build, records[build], flush=True)


if __name__ == "__main__":
    main()
