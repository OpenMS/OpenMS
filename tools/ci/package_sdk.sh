#!/usr/bin/env bash
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
#
# Package the OpenMS developer SDK from a configured and built OpenMS tree and
# check it with the example consumer project src/tests/external.
#
# Usage: package_sdk.sh <build dir> <version> <output dir>
#
# The SDK is the development installation of the core and CLI layers of the
# package (see cmake/install_macros.cmake): their libraries and headers, the CMake
# package files, the shared data and the runtime dependencies of those libraries
# (from the component Dependencies, which exists when the tree was configured
# with a PACKAGE_TYPE, as for the installers, reduced by sdk_prune_dependencies.cmake
# to what the SDK loads). The GUI layer is left out, so the SDK needs no Qt. It is written
# as <output dir>/OpenMS-SDK-<version>-<platform>.tar.gz (.zip on Windows) with a
# single top-level directory of the same name.
#
# Before the archive is accepted, it is extracted to a different location and
# src/tests/external is configured against it through CMAKE_PREFIX_PATH (the way
# the SDK README tells users to), built and its tests are run. The consumer gets
# only the toolchain of the OpenMS build (compiler, vcpkg), which is how it
# obtains the Boost headers when the OpenMS package requires them.
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo >&2 "Usage: $0 <build dir> <version> <output dir>"
  exit 2
fi

build_dir=$(cd "$1" && pwd)
version=$2
mkdir -p "$3"
out_dir=$(cd "$3" && pwd)
source_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cache="$build_dir/CMakeCache.txt"

if [[ ! -f "$cache" ]]; then
  echo >&2 "ERROR: $build_dir is not a configured CMake build tree"
  exit 1
fi

# value of a CMake cache entry (empty if absent)
cache_var() {
  sed -n "s/^$1:[A-Z]*=//p" "$cache" | head -n 1
}

# paths handed to CMake: forward slashes on Windows
cmake_path() {
  if command -v cygpath >/dev/null 2>&1; then
    cygpath -m "$1"
  else
    printf '%s\n' "$1"
  fi
}

case "$(uname -s)" in
  Linux)                platform="Linux-$(uname -m)"; archive_ext=tar.gz ;;
  Darwin)               platform="macOS-$(uname -m)"; archive_ext=tar.gz ;;
  MINGW*|MSYS*|CYGWIN*) platform="Windows-x64";       archive_ext=zip ;;
  *) echo >&2 "ERROR: unsupported system $(uname -s)"; exit 1 ;;
esac

sdk_name="OpenMS-SDK-${version}-${platform}"
work_dir="$out_dir/.sdk-work"
stage_root="$work_dir/stage"
stage="$stage_root/$sdk_name"
archive="$out_dir/$sdk_name.$archive_ext"

rm -rf "$work_dir" "$archive"
mkdir -p "$stage"

#------------------------------------------------------------------------------
# Install the development components into the staging prefix

components=(
  library library_cli
  OpenMS_headers OpenMS_CLI_headers OpenSwathAlgo_headers thirdparty_headers
  cmake cmake_cli
  share
)
# the runtime dependencies of all installed targets, pruned below to those of the
# SDK libraries; their macOS install-name fix-up is redone below, as the
# component's own install code does not work in this layout
components+=(Dependencies)

build_type=$(cache_var CMAKE_BUILD_TYPE)
for component in "${components[@]}"; do
  echo "--- installing component $component"
  cmake --install "$(cmake_path "$build_dir")" --config "${build_type:-Release}" \
        --prefix "$(cmake_path "$stage")" --component "$component" --strip
done

cmake_dir=$(cache_var INSTALL_CMAKE_DIR)
lib_dir=$(cache_var INSTALL_LIB_DIR)
include_dir=$(cache_var INSTALL_INCLUDE_DIR)
for required in \
    "$cmake_dir/OpenMSConfig.cmake" \
    "$cmake_dir/OpenMSConfigVersion.cmake" \
    "$include_dir/OpenMS" \
    "$lib_dir"; do
  if [[ ! -e "$stage/$required" ]]; then
    echo >&2 "ERROR: the SDK misses $required"
    exit 1
  fi
done

# Linux: the SDK has to find its bundled dependencies wherever it is extracted.
# The OpenMS libraries are installed with $ORIGIN/../lib/, but the bundled
# third-party libraries keep whatever RPATH their own build gave them, and with
# RUNPATH the dependencies of a library are looked up through that library's own
# entry only. So every library in the SDK gets one that points at its own
# directory, rather than relying on each dependency provider to have set one.
if [[ "$(uname -s)" == Linux ]]; then
  if ! command -v patchelf >/dev/null 2>&1; then
    echo >&2 "ERROR: patchelf is required to package the SDK on Linux"
    exit 1
  fi
  # shellcheck disable=SC2016 # $ORIGIN is for the dynamic loader, not the shell
  find "$stage/$lib_dir" -maxdepth 1 -type f -name '*.so*' -print0 |
    xargs -0 -r -n 1 patchelf --set-rpath '$ORIGIN'
fi

# macOS: the install code of the Dependencies component (package_mac_productbuild.cmake)
# is written for the pkg layout: its fix_dependencies.rb call names the bin/ directory
# of the application bundle, which does not exist in this prefix, and the script
# aborts on it before it gets to the libraries. So the bundled dependencies keep the
# install names of the build machine. Run the call the library components use
# (-l only) once more over the complete lib/ directory, with @loader_path as the
# only RPATH so the libraries find each other wherever the SDK is extracted (as
# $ORIGIN does on Linux), then sign again: rewriting install names invalidates
# signatures, and arm64 refuses to load such libraries.
if [[ "$(uname -s)" == Darwin ]]; then
  ruby "$source_dir/cmake/MacOSX/fix_dependencies.rb" -l "$stage/$lib_dir" -e @rpath/ -n -c -r @loader_path
  signing_identity=$(cache_var SIGNING_IDENTITY)
  if [[ -n "$signing_identity" && "$signing_identity" != "-" ]]; then
    sign_args=(--force --options runtime --timestamp --sign "$signing_identity")
  else
    sign_args=(--force --sign -)
  fi
  while IFS= read -r -d '' file; do
    if file -b "$file" | grep -q '^Mach-O'; then
      codesign "${sign_args[@]}" "$file"
    fi
  done < <(find "$stage/$lib_dir" -type f -print0)
fi

# Keep only the bundled libraries the SDK libraries load (no Qt, nothing only the
# GUI library or the applications need). The loader environment of this job must
# not decide where a dependency resolves to.
case "$(uname -s)" in
  Linux)  lib_pattern='lib%s.so' ;;
  Darwin) lib_pattern='lib%s.dylib' ;;
  *)      lib_pattern='%s.dll' ;;
esac
roots=()
for library in OpenMS OpenSwathAlgo OpenMS_CLI; do
  # shellcheck disable=SC2059 # the pattern is chosen above
  root="$stage/$lib_dir/$(printf "$lib_pattern" "$library")"
  if [[ ! -e "$root" ]]; then
    echo >&2 "ERROR: the SDK misses $root"
    exit 1
  fi
  roots+=("$(cmake_path "$root")")
done
echo "--- pruning the bundled dependencies"
env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH -u DYLD_FALLBACK_LIBRARY_PATH \
  cmake "-DLIB_DIR=$(cmake_path "$stage/$lib_dir")" "-DROOTS=$(IFS=';'; echo "${roots[*]}")" \
        -P "$(cmake_path "$source_dir/tools/ci/sdk_prune_dependencies.cmake")"

# macOS: the oldest macOS the SDK runs on is the deployment target (minos) its
# libraries were built for. The README states the one of libOpenMS. A bundled
# dependency built for a newer macOS would raise it without anyone noticing, so
# it stops the packaging instead.
macos_minimum=""
if [[ "$(uname -s)" == Darwin ]]; then
  # minos of a Mach-O file: LC_BUILD_VERSION, or LC_VERSION_MIN_MACOSX in older
  # files (reads all of otool's output: an early exit would fail the pipeline)
  macho_minos() {
    otool -l "$1" | awk '
      $1 == "cmd" { cmd = $2 }
      minos == "" && cmd == "LC_BUILD_VERSION" && $1 == "minos" { minos = $2 }
      minos == "" && cmd == "LC_VERSION_MIN_MACOSX" && $1 == "version" { minos = $2 }
      END { print minos }'
  }
  # whether the dotted version $1 is newer than $2
  version_newer() {
    awk -v a="$1" -v b="$2" 'BEGIN {
      n = split(a, x, "."); m = split(b, y, "."); if (m > n) n = m
      for (i = 1; i <= n; ++i) if (x[i] + 0 != y[i] + 0) exit !(x[i] + 0 > y[i] + 0)
      exit 1
    }'
  }
  openms_minimum=$(macho_minos "$stage/$lib_dir/libOpenMS.dylib")
  if [[ -z "$openms_minimum" ]]; then
    echo >&2 "ERROR: libOpenMS.dylib records no minimum macOS version"
    exit 1
  fi
  newer=()
  while IFS= read -r -d '' file; do
    if file -b "$file" | grep -q '^Mach-O'; then
      minimum=$(macho_minos "$file")
      if [[ -n "$minimum" ]] && version_newer "$minimum" "$openms_minimum"; then
        newer+=("${file#"$stage/"}: macOS $minimum")
      fi
    fi
  done < <(find "$stage/$lib_dir" -type f -print0)
  if [[ ${#newer[@]} -gt 0 ]]; then
    echo >&2 "ERROR: libOpenMS is built for macOS $openms_minimum, but these bundled libraries need a newer macOS:"
    printf >&2 '  %s\n' "${newer[@]}"
    exit 1
  fi
  macos_minimum=", deployment target macOS $openms_minimum or newer"
fi

# The README names what a consumer needs, read off the installed package file:
# while the public headers include Boost, OpenMSConfig.cmake records the Boost
# version the build used in _openms_boost_version (falling back to 1.81.0 when the
# build recorded none) and requires it from the consumer; once Boost is a private
# dependency of the shared library, the file has no such entry and a consumer
# needs no Boost development files at all.
config_file="$stage/$cmake_dir/OpenMSConfig.cmake"
if grep -q '^set(_openms_boost_version ' "$config_file"; then
  boost_version=$(sed -n 's/^set(_openms_boost_version "\(.*\)")$/\1/p' "$config_file")
  boost_requirement="Boost headers, version ${boost_version:-1.81.0} or newer (the version this SDK was built
    against; find_package(OpenMS) requires at least that): public OpenMS headers
    include Boost. No compiled Boost library is needed."
else
  boost_requirement="No Boost: the OpenMS headers do not include it, and find_package(OpenMS)
    does not look for it."
fi
awk -v requirement="$boost_requirement" -v major_minor="${version%.*}" -v macos_minimum="$macos_minimum" '
  index($0, "@BOOST_REQUIREMENT@") { sub(/@BOOST_REQUIREMENT@/, requirement) }
  { gsub(/@OPENMS_VERSION_MAJOR_MINOR@/, major_minor); gsub(/@MACOS_MINIMUM@/, macos_minimum); print }
' "$source_dir/cmake/OpenMSSDKReadme.txt" > "$stage/README.txt"
cp "$source_dir/License.txt" "$stage/License.txt"

#------------------------------------------------------------------------------
# Create the archive

echo "--- creating $archive"
if [[ "$archive_ext" == zip ]]; then
  (cd "$stage_root" && cmake -E tar cf "$(cmake_path "$archive")" --format=zip "$sdk_name")
else
  (cd "$stage_root" && cmake -E tar czf "$archive" "$sdk_name")
fi

#------------------------------------------------------------------------------
# Check the archive: extract it elsewhere and build the example consumer against it

check_root="$work_dir/check"
consumer_build="$work_dir/consumer"
mkdir -p "$check_root"
(cd "$check_root" && cmake -E tar xf "$(cmake_path "$archive")")
sdk_prefix="$check_root/$sdk_name"

prefix_path="$(cmake_path "$sdk_prefix")"
build_prefix_path=$(cache_var CMAKE_PREFIX_PATH)
if [[ -n "$build_prefix_path" ]]; then
  prefix_path="$prefix_path;$build_prefix_path"
fi

configure_args=(
  -S "$(cmake_path "$source_dir/src/tests/external")"
  -B "$(cmake_path "$consumer_build")"
  -G "$(cache_var CMAKE_GENERATOR)"
  "-DCMAKE_BUILD_TYPE=${build_type:-Release}"
  "-DCMAKE_PREFIX_PATH=$prefix_path"
)
# the toolchain of the OpenMS build, as the installed-consumer tests forward it
# (src/tests/CMakeLists.txt)
for var in CMAKE_TOOLCHAIN_FILE VCPKG_INSTALLED_DIR VCPKG_TARGET_TRIPLET VCPKG_HOST_TRIPLET \
           CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_MAKE_PROGRAM CMAKE_MSVC_RUNTIME_LIBRARY \
           CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT; do
  value=$(cache_var "$var")
  if [[ -n "$value" ]]; then
    configure_args+=("-D$var=$value")
  fi
done

echo "--- building src/tests/external against the SDK"
cmake "${configure_args[@]}"

# the consumer has to have picked up the SDK, not some other OpenMS
found_dir=$(sed -n 's/^OpenMS_DIR:PATH=//p' "$consumer_build/CMakeCache.txt")
case "$found_dir" in
  "$(cmake_path "$sdk_prefix")"/*) ;;
  *) echo >&2 "ERROR: the consumer found OpenMS at '$found_dir', not in the SDK"; exit 1 ;;
esac

cmake --build "$(cmake_path "$consumer_build")" --config "${build_type:-Release}"
# as the README tells users to run their programs
OPENMS_DATA_PATH="$(cmake_path "$sdk_prefix/$(cache_var INSTALL_SHARE_DIR)")"
export OPENMS_DATA_PATH
ctest --test-dir "$(cmake_path "$consumer_build")" -C "${build_type:-Release}" --output-on-failure --no-tests=error

rm -rf "$work_dir"
echo "SDK archive: $archive"
