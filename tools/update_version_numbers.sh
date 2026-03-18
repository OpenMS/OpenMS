#!/usr/bin/env bash

################################################################################
set -eu
set -o pipefail

################################################################################
usage() {
  cat <<EOF
Usage: $(basename "$0") [options] major minor patch

  -h      This message

Update the OpenMS version number.

EOF
}

################################################################################
while getopts "h" o; do
  case "${o}" in
  h)
    usage
    exit
    ;;

  *)
    exit 1
    ;;
  esac
done

shift $((OPTIND - 1))

if [ $# -ne 3 ]; then
  echo >&2 "ERROR: please provide three numbers (use --help for more info)"
  exit 1
fi

################################################################################
package_version_major=$1
package_version_minor=$2
package_version_patch=$3
package_version="$1.$2.$3"

echo "Setting version $package_version"

################################################################################
# update main cmakelist
sed -i -e "s#.*set(OPENMS_PACKAGE_VERSION_MAJOR.*#set(OPENMS_PACKAGE_VERSION_MAJOR \"$package_version_major\")#" CMakeLists.txt
sed -i -e "s#.*set(OPENMS_PACKAGE_VERSION_MINOR.*#set(OPENMS_PACKAGE_VERSION_MINOR \"$package_version_minor\")#" CMakeLists.txt
sed -i -e "s#.*set(OPENMS_PACKAGE_VERSION_PATCH.*#set(OPENMS_PACKAGE_VERSION_PATCH \"$package_version_patch\")#" CMakeLists.txt

# update version info test
sed -i -e "s#detail.version_major =.*#detail.version_major = $package_version_major;#" ./src/tests/class_tests/openms/source/VersionInfo_test.cpp
sed -i -e "s#detail.version_minor =.*#detail.version_minor = $package_version_minor;#" ./src/tests/class_tests/openms/source/VersionInfo_test.cpp
sed -i -e "s#detail.version_patch =.*#detail.version_patch = $package_version_patch;#" ./src/tests/class_tests/openms/source/VersionInfo_test.cpp

# update vcpkg.json
sed -i -e "s/\"version-string\": \".*\"/\"version-string\": \"$package_version\"/" vcpkg.json

# update test write ini out:
sed -i -e "s#<ITEM name=\"version\" value=\".*\" type=\"string\"#<ITEM name=\"version\" value=\"$package_version\" type=\"string\"#g" ./src/tests/topp/WRITE_INI_OUT.ini
# update INIUpdater version
sed -i -e "s#<ITEM name=\"version\" value=\".*\" type=\"string\"#<ITEM name=\"version\" value=\"$package_version\" type=\"string\"#g" ./src/tests/topp/INIUpdater_1_noupdate.toppas

# update INIs in tests topp:
find ./src/tests/topp/ \
  -type f \
  -name '*.ini' \
  -exec grep -q "<ITEM name=\"version\" value=\".*\" type=\"string\"" {} \; \
  -exec sed -i -e "s#<ITEM name=\"version\" value=\".*\" type=\"string\"#<ITEM name=\"version\" value=\"$package_version\" type=\"string\"#g" {} \;

# Update Changelog
#
# Create a banner with line endings that sed wants (add a backslash
# in front of each newline except the last.
section_header=$(
  sed -Ee '$!s/$/\\/' <<<"
------------------------------------------------------------------------------------------
----                                OpenMS ${package_version}     (under development)              ----
------------------------------------------------------------------------------------------

 " # NOTE: this line should start with a single blank space!
)

# Get the line number of the first line starting with more than 10 "-"
line_num=$(
  grep \
    --extended-regexp \
    --line-number \
    --max-count=1 "^-{10,}" CHANGELOG |
    cut -d: -f1
)

sed -i \
  -e "${line_num}i $section_header" \
  CHANGELOG

################################################################################
echo "done."
