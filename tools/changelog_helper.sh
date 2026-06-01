#!/bin/bash
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Johannes Veit $
# $Authors: Johannes Veit $
# --------------------------------------------------------------------------

############################################################################
#
# Compare the binary folders of two OpenMS releases (old / new) and output:
#
# - Discontinued tools
# - New tools
# - Parameter changes from generated INI files
#
# Output is in Markdown format.
#
# Worked for me on MacOS comparing release 2.0 (built from source) vs.
# release 1.11.1 (installed using binary installer). Might not run
# out-of-the-box on linux because of differences in sed syntax etc.
#
############################################################################

# check command line arguments
if [[ $# != 2 ]]
then 
    echo "Usage: $0 <old-release-bin-dir> <new-release-bin-dir>"
    exit 1
fi

# directories containing TOPP binaries
BIN_DIR_OLD=$1
BIN_DIR_NEW=$2

# tmp files
SYSTEM_TMP_DIR=/tmp
TMP_DIR=${SYSTEM_TMP_DIR}/OpenMS_changelog_helper
mkdir ${TMP_DIR}
mkdir ${TMP_DIR}/inis
mkdir ${TMP_DIR}/inis/old
mkdir ${TMP_DIR}/inis/new
TMP_FILE_OLD=${TMP_DIR}/tool_list_old.txt
TMP_FILE_NEW=${TMP_DIR}/tool_list_new.txt
TMP_FILE_COMM=${TMP_DIR}/common_tools.txt

# store relevant tool names in tmp files (remove e.g., GUI tools)
ls -la ${BIN_DIR_OLD}/ \
    | awk '{print $9}' \
    | sort \
    | grep -v -e "Tutorial\|TOPPAS\|TOPPView\|INIFileEditor\|SEARCHENGINES\|OpenMSInfo\|GenericWrapper\|SwathWizard\|FLASHDeconvWizard\|Testing" \
    | grep -v -e "\.$" \
    | grep -v -e "^$" \
    > ${TMP_FILE_OLD}

ls -la ${BIN_DIR_NEW}/ \
    | awk '{print $9}' \
    | sort \
    | grep -v -e "Tutorial\|TOPPAS\|TOPPView\|INIFileEditor\|SEARCHENGINES\|OpenMSInfo\|GenericWrapper\|SwathWizard\|FLASHDeconvWizard\|Testing" \
    | grep -v -e "\.$" \
    | grep -v -e "^$"  \
    > ${TMP_FILE_NEW}

# find removed tools and new tools
for s in ADDED REMOVED
do
    BIN_DIR=${BIN_DIR_NEW}
    GREP_CHAR=">"
    if [[ $s == "REMOVED" ]]
    then
        BIN_DIR=${BIN_DIR_OLD}
        GREP_CHAR="<"
    fi
    echo
    echo "## $s"
    echo
    diff ${TMP_FILE_OLD} ${TMP_FILE_NEW} \
        | grep -e "^${GREP_CHAR}" \
        | sort \
        | while read i
        do
            TOOL_NAME=$(echo $i | sed -E "s/^. //")
            TOOL_DESCR=$(LD_LIBRARY_PATH=${BIN_DIR}/../lib:${LD_LIBRARY_PATH} ${BIN_DIR}/${TOOL_NAME} --help 2>&1 | grep " -- " | head -n 1 | sed -E 's/.* -- (.*)$/\1/' | sed -E 's/\.$//')
            if [[ ${TOOL_DESCR} != "" ]]
            then
                echo "- ${TOOL_NAME} -- ${TOOL_DESCR}"
            else
                echo "- ${TOOL_NAME}"
            fi


        done
    echo
done

# store names of tools present in both old and new release in tmp file
comm -12 ${TMP_FILE_OLD} ${TMP_FILE_NEW} > ${TMP_FILE_COMM}

# print changed parameters as markdown table
echo
echo "## CHANGED PARAMETERS"
echo
echo "| Tool name | Added/removed | Parameter name | Type | Default value | Restrictions | Supported formats |"
echo "|-----------|:-------------:|----------------|------|---------------|--------------|-------------------|"

# write ini files for old and new tools, modify them on the fly:
#
# - remove stuff where we're not interested in changes
# - replace parameter names with:this:notation for nested parameters
#
# (=> result is not a valid INI file anymore)
for INI_SUB_DIR in old new
do
    BIN_DIR=${BIN_DIR_NEW}
    if [[ ${INI_SUB_DIR} == "old" ]]
    then
        BIN_DIR=${BIN_DIR_OLD}
    fi

    cat ${TMP_FILE_COMM} | while read t
    do
        # reset subsection prefix
        p=
        # generate pseudo ini file
        LD_LIBRARY_PATH=${BIN_DIR}/../lib:${LD_LIBRARY_PATH} ${BIN_DIR}/$t -write_ini - \
            | grep -v "name=\"version\"" \
            | sed -E 's/description="[^"]*"|required="[^"]*"|advanced="[^"]*"//g' \
            | sed -E 's/restrictions="[^"]*pyrophospho[^"]*"/ restrictions="..."/' \
            | sed -E 's/tags="is_executable"//g' \
            | while read l
            do
                # new NODE -> append subsection to prefix p
                echo $l | grep "<NODE"  &> /dev/null && p=$p$(echo $l | sed -E 's/.*name="([^"]+)".*/\1/'): && continue
                # NODE finished -> remove last subsection from prefix p
                echo $l | grep "</NODE" &> /dev/null && p=$(echo $p | sed -E 's/[^:]+:$//') && continue
                # otherwise, substitute name -> name:with:subsections
                echo $l | sed -E 's/name="([^"]*)"/name="'$p'\1"/'
            done \
            > ${TMP_DIR}/inis/${INI_SUB_DIR}/$t.pseudo.ini 2> /dev/null
    done
done

# sort pseudo ini file
cat ${TMP_FILE_COMM} | while read t
do
    sort ${TMP_DIR}/inis/old/$t.pseudo.ini -o ${TMP_DIR}/inis/old/$t.pseudo.ini
    sort ${TMP_DIR}/inis/new/$t.pseudo.ini -o ${TMP_DIR}/inis/new/$t.pseudo.ini
done \

# compute diffs of pseudo ini files and output markdown table of changed parameters
cat ${TMP_FILE_COMM} | while read t
do
    diff -d ${TMP_DIR}/inis/old/$t.pseudo.ini ${TMP_DIR}/inis/new/$t.pseudo.ini
done \
    | grep "<ITEM" \
    | sed -E 's/^<[[:space:]]+/- /' \
    | sed -E 's/^>[[:space:]]+/+ /' \
    | perl -0777 -pe 's/\n+/\n/g' \
    | while read l
    do
        T_NAME=$(echo $l | sed -E 's/.*name="([^":]*):.*/\1/')
        P_ADD_REM=$(echo $l | sed -E 's/^(.).*/\1/')
        P_NAME=$(echo $l | sed -E 's/.*name="[^":]+:.:([^"]*)".*/\1/')
        P_TYPE=$(echo $l | grep "type=" | sed -E 's/.*type="([^"]*)".*/\1/')
        P_VALUE=$(echo $l | grep "value=" | sed -E 's/.*value="([^"]*)".*/\1/')
        P_RESTRICTIONS=$(echo $l | grep "restrictions=" | sed -E 's/.*restrictions="([^"]*)".*/\1/')
        P_FORMATS=$(echo $l | grep "supported_formats=" | sed -E 's/.*supported_formats="([^"]*)".*/\1/')
        echo "| ${T_NAME} | ${P_ADD_REM} | ${P_NAME} | ${P_TYPE} | ${P_VALUE} | ${P_RESTRICTIONS} | ${P_FORMATS} |"
    done | sort

#cleanup
rm -rf ${TMP_DIR}

