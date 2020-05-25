#!/bin/bash

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

# This script finds all the required qt libs for the OpenMS installation

# we assume to have 3 arguments
# 1. an executable linking against all shipped libraries
# 2. the path where the OpenMS libs where build
# 3. the path were the external libs should be placed

if [ ! $# == 3 ]; then
  echo "Usage: $0 /path/to/executables /path/to/libs /target/path"
  exit
fi

ldd_targets=$(find $1 $2 -type f)
target_path=$3

globalSOs=()

for exe in $ldd_targets; do
  currentSOs=$(ldd $exe | grep "libQt" | sed 's/.* => \(.*\) .*/\1/g')

  for so in ${currentSOs}; do
    #echo "cp ${so} ${target_path}"
    #globalSOs$=($(printf "%s\n" "${globalSOs[@]}" | sort -u))
    globalSOs+=($so)
  done
done

# make list unique
uniqSOs=$(echo "${globalSOs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')

for so in ${uniqSOs[@]}; do
  cp ${so} ${target_path}
done
