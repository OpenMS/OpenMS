#!/usr/bin/env bash

################################################################################
set -eu
set -o pipefail

################################################################################
vars_to_cache=(
  "ADDRESS_SANITIZER"
  "CC"
  "CXX"
  "CMAKE_CXX_COMPILER"
  "CFLAGS"
  "CXXFLAGS"
  "CMAKE_PREFIX_PATH"
  "CMAKE_BUILD_TYPE"
  "CMAKE_GENERATOR_PLATFORM"
  "Boost_DEBUG"
  "BOOST_USE_STATIC"
  "OPENMS_CONTRIB_LIBS"
  "ENABLE_CLASS_TESTING"
  "ENABLE_GCC_WERROR"
  "ENABLE_STYLE_TESTING"
  "ENABLE_TOPP_TESTING"
  "ENABLE_PIPELINE_TESTING"
  "ENABLE_DOCS"
  "ENABLE_CWL_GENERATION"
  "ENABLE_TUTORIALS"
  "ENABLE_UPDATE_CHECK"
  "MT_ENABLE_OPENMP"
  "NO_DEPENDENCIES"
  "SEARCH_ENGINES_DIRECTORY"
  "SIGNING_IDENTITY"
  "PACKAGE_TYPE"
  "PYOPENMS"
  "PY_MEMLEAK_DISABLE"
  "PY_NO_OPTIMIZATION"
  "PY_NO_OUTPUT"
  "PY_NUM_MODULES"
  "PY_NUM_THREADS"
  "Python_ROOT_DIR"
  "Python_FIND_STRATEGY"
  "WITH_GUI"
  "WITH_PARQUET"
  "WITH_THERMORAWFILEPARSER_TEST"
  "COMPILE_PXDS"
  "USE_EXTERNAL_JSON"
)

################################################################################
option_verbose=0

################################################################################
function usage() {
  cat <<EOF
Usage: $(basename "$0") [options] file

  -h      This message
  -v      Verbose output

Write environment variables to file.
EOF
}

################################################################################
function write_file() {
  local file=$1

  mkdir -p "$(dirname "$file")"

  if [ -e "$file" ]; then
    rm -f "$file"
  fi

  for var in "${vars_to_cache[@]}"; do
    # This next line is a little tricky.  It says to evaluate $var and
    # use its value as the name of another variable to evaluate
    # (indirection).  The ":-" bit says to set the variable to an
    # empty string if it is not already defined.
    val="${!var:-}"

    if [ -n "${val}" ]; then
      if [ "$option_verbose" -eq 1 ]; then
        echo "Found $var with value $val"
      fi

      printf '%s:STRING=%s\n' "$var" "${val}" >>"$file"
    fi
  done
}

################################################################################
function main() {
  while getopts "hv" o; do
    case "${o}" in
    h)
      usage
      exit
      ;;

    v)
      option_verbose=1
      ;;

    *)
      exit 1
      ;;
    esac
  done

  shift $((OPTIND - 1))

  if [ $# -ne 1 ]; then
    echo >&2 "ERROR: missing file argument, see -h"
    exit 1
  fi

  write_file "$1"
}

################################################################################
main "$@"
