#!/bin/bash
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
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

# directories containing TOPP binaries
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
    | grep -v -e "Tutorial\|TOPPAS\|TOPPView\|INIFileEditor\|SEARCHENGINES\|OpenMSInfo\|GenericWrapper\|Testing\|SwathWizard" \
    | grep -v -e "\.$" \
    | grep -v -e "^$"  \
    | while read i
    do
      echo -e $i'\t'$(${BIN_DIR}/${i} --help 2>&1 | grep "Version" | sed -E "s/Version: //" | sed -E 's/\s.*$/ /')
   done

