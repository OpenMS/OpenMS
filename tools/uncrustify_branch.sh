#!/usr/bin/env sh

# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2020.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

# script to apply uncrustify to all files modified in the
# current branch

umask 002
SCRIPT_PATH=`dirname $0`

usage()
{
  cat << EOF
usage: $0 options

This script runs uncrustify on all files modified in the current branch.

OPTIONS:
   -h      Show this message
   -r      Refence branch [default: develop]
   -f      Perform uncrustify even if the current HEAD is dirty (files where modified).
   -n      Dry-run. Only show what would be uncrustify(ed).
   -b      Create separate branch.
   -v      Verbose
EOF
}

print_if_verbose()
{
  if [ $VERBOSE -ne 0 ]; then
    echo "-> $1"
  fi
}

REF_BRANCH=
VERBOSE=0
FORCE=0
DRY=0
CREATE_BRANCH=0

while getopts “hr:vfnb” OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    r)
      REF_BRANCH=$OPTARG
      ;;
    f)
      FORCE=1
      ;;
    v)
      VERBOSE=1
      ;;
    n)
      DRY=1
      ;;
    b)
      CREATE_BRANCH=1
      ;;
    ?)
      usage
      exit
      ;;
  esac
done

# set ref branch to develop if none was given
if [[ -z $REF_BRANCH ]]
then
  REF_BRANCH=develop
fi

# check if uncrustify is installed
UNCRUCSTIFY=$(which uncrustify)
if [ $? -ne 0 ]; then
  echo "Missing uncrustify on this system! Exit!"
  exit 1
fi

# go to OpenMS dir ${SCRIPT_PATH}/..
OpenMS_DIR=${SCRIPT_PATH}/..
pushd ${OpenMS_DIR} > /dev/null

if [ ! -e `pwd`/tools/uncrustify.cfg ]; then
  echo "Missing uncrustify.cfg! Is this really an OpenMS clone? Exit!"
  exit 1
fi

# check if branch is dirty
git diff-index --quiet HEAD
if [ $? -ne 0 ] && [ $FORCE -ne 1 ]; then
  echo "Your branch contains modified files. Please commit them before applying this script to avoid mixing automated- and user-changes."
  echo "If you know what you are doing you can turn this check of with -f".
  exit 1
fi

# get the current branch name
branch_name=$(git rev-parse --abbrev-ref HEAD)

print_if_verbose "Applying to branch ${branch_name}"

# get modified files .. including only *.cpp and *.h, but no tests
modified_files=$(git diff --name-only ${branch_name} $(git merge-base ${branch_name} ${REF_BRANCH}) | grep -e '.cpp$' -e '.h$' | grep -v '_test.cpp$' | grep -v 'thirdparty')

# create uncrustified branch 
if [ $CREATE_BRANCH == 1 ]; then
  uncrustify_branch_name=${branch_name}_uncrustify
  echo "Creating and switching to uncrustify branch ${uncrustify_branch_name}"
  echo "Run git commit -a to commit uncrustified files afterwards!"

  if [ $DRY -ne 1 ]; then
    git checkout -b ${uncrustify_branch_name} ${branch_name}
  else
    echo "git checkout -b ${uncrustify_branch_name} ${branch_name}"
  fi
fi

for f in $modified_files; do
  print_if_verbose "Uncrustify file: $f"
  if [ $DRY -ne 1 ]; then
    ${UNCRUCSTIFY} -lcpp -q -c tools/uncrustify.cfg --no-backup $f
  else
    echo "${UNCRUCSTIFY} -lcpp -q -c tools/uncrustify.cfg --no-backup $f"
  fi
done

popd > /dev/null
