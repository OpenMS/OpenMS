#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
"""Fail if a macOS wheel needs a newer macOS than its platform tag promises.

pip selects a wheel by its tag alone, so a wheel whose bundled dependencies target a
newer macOS installs on an older one and then dies in dyld.

Usage: check_macos_wheel_target.py wheelhouse/*.whl
"""
import re
import struct
import sys
import zipfile

MH_MAGIC_64, MH_CIGAM_64 = 0xFEEDFACF, 0xCFFAEDFE
FAT_MAGIC, FAT_CIGAM = 0xCAFEBABE, 0xBEBAFECA
LC_VERSION_MIN_MACOSX, LC_BUILD_VERSION = 0x24, 0x32


def decode(version):
    return (version >> 16, (version >> 8) & 0xFF, version & 0xFF)


def min_os(blob):
    """Highest minimum-macOS requirement in one Mach-O (or fat) image, or None."""
    if len(blob) < 8:
        return None
    magic = struct.unpack(">I", blob[:4])[0]
    if magic in (FAT_MAGIC, FAT_CIGAM):
        count = struct.unpack(">I", blob[4:8])[0]
        found = []
        for index in range(count):
            start = 8 + index * 20
            if start + 20 > len(blob):
                break
            offset, size = struct.unpack(">II", blob[start + 8:start + 16])
            slice_version = min_os(blob[offset:offset + size])
            if slice_version:
                found.append(slice_version)
        return max(found) if found else None

    little = struct.unpack("<I", blob[:4])[0] == MH_MAGIC_64
    big = struct.unpack(">I", blob[:4])[0] == MH_MAGIC_64
    if not (little or big):
        return None
    endian = "<" if little else ">"
    ncmds = struct.unpack(endian + "I", blob[16:20])[0]
    offset, best = 32, None
    for _ in range(ncmds):
        if offset + 8 > len(blob):
            break
        cmd, cmdsize = struct.unpack(endian + "II", blob[offset:offset + 8])
        if cmdsize < 8:
            break
        if cmd == LC_BUILD_VERSION and offset + 16 <= len(blob):
            best = max(best or (0, 0, 0), decode(struct.unpack(endian + "I", blob[offset + 12:offset + 16])[0]))
        elif cmd == LC_VERSION_MIN_MACOSX and offset + 12 <= len(blob):
            best = max(best or (0, 0, 0), decode(struct.unpack(endian + "I", blob[offset + 8:offset + 12])[0]))
        offset += cmdsize
    return best


def check(path):
    tag = re.search(r"-macosx_(\d+)_(\d+)_", path)
    if not tag:
        print(f"  skipped (not a macOS wheel): {path}")
        return True
    promised = (int(tag.group(1)), int(tag.group(2)), 0)
    offenders = []
    with zipfile.ZipFile(path) as wheel:
        for item in wheel.namelist():
            if item.endswith("/"):
                continue
            with wheel.open(item) as handle:
                head = handle.read(4)
                if len(head) < 4 or struct.unpack(">I", head)[0] not in (
                        MH_MAGIC_64, MH_CIGAM_64, FAT_MAGIC, FAT_CIGAM):
                    continue
            required = min_os(wheel.read(item))
            if required and required > promised:
                offenders.append((item, required))
    name = path.rsplit("/", 1)[-1]
    if offenders:
        print(f"  FAIL {name}: tag promises macOS "
              f"{promised[0]}.{promised[1]} but {len(offenders)} binaries need more")
        for item, required in sorted(offenders, key=lambda entry: entry[1], reverse=True)[:10]:
            print(f"        {required[0]}.{required[1]} {item}")
        return False
    print(f"  ok   {name}: every binary runs on macOS {promised[0]}.{promised[1]}")
    return True


def main(paths):
    if not paths:
        print("No wheels given", file=sys.stderr)
        return 2
    print(f"Checking {len(paths)} wheel(s) against their macOS platform tag")
    if all([check(path) for path in paths]):
        return 0
    print("\nA wheel would install on a macOS it cannot run on. Either raise\n"
          "MACOSX_DEPLOYMENT_TARGET to match the bundled dependencies, or build the\n"
          "dependencies for the lower target (a runner image matching it).", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
