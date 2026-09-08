#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $

"""Install an existing CI build and compile/run the external CMake consumer."""

import argparse
import os
from pathlib import Path
import subprocess
import tempfile


def read_cache(path):
    """Read CMake cache values, preserving spaces and embedded equals signs."""
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("#", "//")) or "=" not in line:
            continue
        key, value = line.split("=", 1)
        name, separator, _ = key.partition(":")
        if separator:
            values[name] = value
    return values


def run(command, env=None):
    print(subprocess.list2cmdline([str(arg) for arg in command]), flush=True)
    subprocess.run(command, check=True, env=env)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_dir", type=Path, help="Configured and already built OpenMS CI directory")
    parser.add_argument("--config", default="Release", help="Configuration for a multi-config build")
    args = parser.parse_args()
    build = args.build_dir.resolve()
    cache = read_cache(build / "CMakeCache.txt")
    source = Path(cache["CMAKE_HOME_DIRECTORY"])
    prefix = Path(cache["CMAKE_INSTALL_PREFIX"])
    config = cache.get("CMAKE_BUILD_TYPE") or args.config
    cmake = cache.get("CMAKE_COMMAND", "cmake")
    ctest = cache.get("CMAKE_CTEST_COMMAND", "ctest")

    # Reuse the completed library build. Avoid application, documentation and
    # redistribution components: this test is a development-package consumer.
    for component in ("library", "OpenMS_headers", "OpenSwathAlgo_headers", "cmake", "share"):
        run([cmake, "--install", str(build), "--config", config, "--component", component])

    package_dir = prefix / cache["INSTALL_CMAKE_DIR"]
    if not (package_dir / "OpenMSConfig.cmake").is_file():
        raise RuntimeError(f"Installed OpenMSConfig.cmake is missing from {package_dir}")

    with tempfile.TemporaryDirectory(prefix="openms-consumer-") as temporary:
        consumer = Path(temporary) / "build"
        command = [
            cmake, "-S", str(source / "src/tests/external"), "-B", str(consumer),
            "-G", cache["CMAKE_GENERATOR"],
            f"-DOpenMS_DIR={package_dir.as_posix()}",
            f"-DOPENMS_EXPECTED_PREFIX={prefix.as_posix()}",
            f"-DCMAKE_BUILD_TYPE={config}",
            f"-DCMAKE_RUNTIME_OUTPUT_DIRECTORY={prefix.as_posix()}/bin",
            f"-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_{config.upper()}={prefix.as_posix()}/bin",
            "-DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF",
            "-DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF",
            "-DCMAKE_DISABLE_FIND_PACKAGE_CURL=ON",
            "-DCMAKE_DISABLE_FIND_PACKAGE_XercesC=ON",
            "-DVCPKG_MANIFEST_MODE=OFF",
        ]
        # Reuse the compiler and dependency installation, without installing a
        # second vcpkg manifest or propagating OpenMS build-tree include paths.
        for key in (
            "CMAKE_TOOLCHAIN_FILE", "VCPKG_INSTALLED_DIR", "VCPKG_TARGET_TRIPLET",
            "VCPKG_HOST_TRIPLET", "CMAKE_PREFIX_PATH", "CMAKE_C_COMPILER",
            "CMAKE_CXX_COMPILER", "CMAKE_MSVC_RUNTIME_LIBRARY", "CMAKE_OSX_ARCHITECTURES",
            "CMAKE_OSX_DEPLOYMENT_TARGET", "CMAKE_OSX_SYSROOT",
        ):
            if cache.get(key):
                command.append(f"-D{key}={cache[key]}")
        for key, flag in (("CMAKE_GENERATOR_PLATFORM", "-A"), ("CMAKE_GENERATOR_TOOLSET", "-T")):
            if cache.get(key):
                command.extend([flag, cache[key]])
        run(command)
        run([cmake, "--build", str(consumer), "--config", config, "--parallel", "2"])

        env = os.environ.copy()
        dependency_prefixes = [Path(p) for p in cache.get("CMAKE_PREFIX_PATH", "").split(";") if p]
        if cache.get("VCPKG_INSTALLED_DIR") and cache.get("VCPKG_TARGET_TRIPLET"):
            dependency_prefixes.append(Path(cache["VCPKG_INSTALLED_DIR"]) / cache["VCPKG_TARGET_TRIPLET"])
        runtime = [prefix / cache["INSTALL_LIB_DIR"], prefix / "bin"]
        for dep in dependency_prefixes:
            if config.lower() == "debug":
                runtime.extend([dep / "debug/bin", dep / "debug/lib"])
            runtime.extend([dep / "bin", dep / "lib"])
        for key in ("PATH", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"):
            env[key] = os.pathsep.join([str(p) for p in runtime] + ([env[key]] if env.get(key) else []))
        env.pop("OPENMS_DATA_PATH", None)

        # File's development fallback otherwise finds the source data directory,
        # hiding missing install rules. Restore it even when the consumer fails.
        source_data = source / "share/OpenMS"
        with tempfile.TemporaryDirectory(prefix="consumer-data-", dir=source_data.parent) as backup:
            saved_data = Path(backup) / "OpenMS"
            source_data.rename(saved_data)
            try:
                run([ctest, "--test-dir", str(consumer), "-C", config,
                     "--output-on-failure", "--no-tests=error"], env=env)
            finally:
                saved_data.rename(source_data)


if __name__ == "__main__":
    main()
