#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Karnajeet Gosavi $

"""Compare object-command durations in saved OFF/ON Ninja logs."""

import argparse
from pathlib import Path


CONVERTED = ("StringUtils.cpp.obj", "DateTime.cpp.obj", "File.cpp.obj")


def read_objects(path):
    """Keep the last compilation of each object; exclude scans and BMI aliases."""
    objects = {}
    for number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 5:
            raise ValueError(f"{path}:{number}: invalid Ninja log row")
        start, end = int(fields[0]), int(fields[1])
        if end < start:
            raise ValueError(f"{path}:{number}: negative command duration")
        output = fields[3].replace("\\", "/")
        if output.endswith((".obj", ".o")):
            objects[output] = (end - start) / 1000
    return objects


def summarize(off, on):
    """Return a Markdown comparison of object durations (in seconds).

    Separate matching outputs from added module objects and converted sources
    from the rest of OpenMS. Raise ValueError if there are no matching outputs.
    Command durations overlap, so their sums do not represent build wall time.
    """
    common = off.keys() & on.keys()
    if not common:
        raise ValueError("The logs have no matching object outputs to compare")
    converted = {p for p in common if "/OpenMS.dir/" in p and p.rsplit("/", 1)[-1] in CONVERTED}
    groups = (
        ("All object commands (including module objects)", off.keys(), on.keys()),
        ("Matching objects", common, common),
        ("Converted sources", converted, converted),
        ("Other matching objects", common - converted, common - converted),
    )
    lines = [
        "| Object commands | OFF count | ON count | OFF sum (s) | ON sum (s) | ON - OFF (s) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for label, off_names, on_names in groups:
        off_sum = sum(off[p] for p in off_names)
        on_sum = sum(on[p] for p in on_names)
        lines.append(f"| {label} | {len(off_names)} | {len(on_names)} | "
                     f"{off_sum:.3f} | {on_sum:.3f} | {on_sum - off_sum:+.3f} |")
    lines.extend([
        "",
        "These sums are elapsed command durations, not CPU time or wall-clock build time; "
        "parallel commands overlap. Scans, module-map generation and linking are excluded. "
        "Only the latest entry per object is counted. Check the build/test outcomes before "
        "interpreting partial logs, and use repeated interleaved clean OFF/ON pairs to assess variability.",
    ])
    return "\n".join(lines)


def main():
    """Read the OFF/ON logs and print a report, or an argparse error on bad input."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("off_log", type=Path)
    parser.add_argument("on_log", type=Path)
    args = parser.parse_args()
    try:
        print(summarize(read_objects(args.off_log), read_objects(args.on_log)))
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
