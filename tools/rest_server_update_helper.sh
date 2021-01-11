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
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

############################################################################
# Scans the binary folder and outputs tool name and version string used in
# the REST update server. Pipe the output in a versions.txt file and 
# replace the the corresponding file on the REST server.
############################################################################

# check command line arguments
if [[ $# != 1 ]]
then 
    echo "Usage: $0 <new-release-bin-dir>"
    exit 1
fi

# directories containing TOPP/UTILS binaries
BIN_DIR=$1

# tmp filesV
SYSTEM_TMP_DIR=/tmp
TMP_DIR=${SYSTEM_TMP_DIR}/OpenMS_REST_update
mkdir ${TMP_DIR}
TMP_FILE_NEW=${TMP_DIR}/tool_list.txt

# store relevant tool names in tmp files
ls -la ${BIN_DIR}/ \
    | awk '{print $9}' \
    | sort \
    | grep -v -e "Tutorial\|TOPPAS\|TOPPView\|INIFileEditor\|SEARCHENGINES\|OpenMSInfo\|GenericWrapper" \
    | grep -v -e "\.$" \
    | grep -v -e "^$"  \
    | while read i
    do
      echo -e $i'\t'$(${BIN_DIR}/${i} --help 2>&1 | grep "Version" | sed -E "s/Version: //" | sed -E 's/\s.*$/ /')
   done

