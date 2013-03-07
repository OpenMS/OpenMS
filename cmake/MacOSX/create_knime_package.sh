#!/bin/bash

# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

## semi-automatic script to prepare the knime package for mac osx given a 
## fully build OpenMS

# check if we are in the correct directory (descriptors etc.) 
if [ ! -e ${PWD}/icons -a ! -e ${PWD}/descriptors ]; then
  echo "Please run the target prepare_knime_package and change to the directory <OpenMS-build>/ctds"
  exit 1
fi

KNIME_PACKAGE=${PWD}/
DESCRIPTORS=${KNIME_PACKAGE}descriptors/
BUILD_DIR=${PWD}/../
PAYLOAD=${KNIME_PACKAGE}payload

SOURCE_PATH=`grep -E "OpenMS_SOURCE_DIR" "${BUILD_DIR}/CMakeCache.txt" | sed 's/.*=//g'`

echo "$SOURCE_PATH"

if [ ! -e ${PAYLOAD} ]; then
  echo "Create missing payload directory"
  mkdir ${PAYLOAD}
fi
  
# script path
SCRIPT_PATH=`dirname "$0"`

echo -n "Copy binaries .. "
cp -r ${BUILD_DIR}bin ${PAYLOAD}
echo "done"

echo -n "Clean up binaries .. "
rm -r ${PAYLOAD}/bin/TOPPView.app
rm -r ${PAYLOAD}/bin/INIFileEditor.app
rm -r ${PAYLOAD}/bin/TOPPAS.app

rm -r ${PAYLOAD}/bin/PhosphoScoring 
rm -r ${PAYLOAD}/bin/OpenMSInfo 
rm -r ${PAYLOAD}/bin/FuzzyDiff 
rm -r ${PAYLOAD}/bin/GenericWrapper

rm -r ${PAYLOAD}/bin/PILISIdentification
rm -r ${PAYLOAD}/bin/PILISModelCV
rm -r ${PAYLOAD}/bin/PILISModelTrainer
rm -r ${PAYLOAD}/bin/PILISSpectraGenerator
rm -r ${PAYLOAD}/bin/SvmTheoreticalSpectrumGeneratorTrainer
rm -r ${PAYLOAD}/bin/InspectAdapter 
rm -r ${PAYLOAD}/bin/MascotAdapter

# Myrimatch is not supported on mac osx
rm -r ${PAYLOAD}/bin/MyriMatchAdapter

# Also invald
rm -r ${PAYLOAD}/bin/OpenSwathMzMLFileCacher

# we do not ship the binary so it will most likely not => remove for now
rm -r ${PAYLOAD}/bin/PepNovoAdapter
echo "done"

echo -n "Copy libraries .. "
cp -r ${BUILD_DIR}lib ${PAYLOAD}
echo "done"

echo -n "Fix libraries and binaries .. "
${SCRIPT_PATH}/fix_dependencies.rb  -l ${PAYLOAD}/lib/ -b ${PAYLOAD}/bin/
echo "done"

echo -n "Fixing ctds .. "

# remove debug output from IDPosteriorErrorProbability
grep -v output_name ${DESCRIPTORS}IDPosteriorErrorProbability.ctd > temp.ctd
mv temp.ctd ${DESCRIPTORS}IDPosteriorErrorProbability.ctd
  
# remove binary path from (omssa|xtandem)Adapter.ctd
grep -v omssa_executable ${DESCRIPTORS}OMSSAAdapter.ctd > omssa.ctd
mv omssa.ctd  ${DESCRIPTORS}OMSSAAdapter.ctd
grep -v xtandem_executable ${DESCRIPTORS}XTandemAdapter.ctd > tandem.ctd
mv tandem.ctd ${DESCRIPTORS}XTandemAdapter.ctd
grep -v default_input_file ${DESCRIPTORS}XTandemAdapter.ctd > tandem.ctd
mv tandem.ctd ${DESCRIPTORS}XTandemAdapter.ctd

echo "done"

echo -n "Export share .. "
svn export ${SOURCE_PATH}/share ${PAYLOAD}/share
echo "done"

echo -n "Copy binaries_mac_64.ini .. "
cp ${SOURCE_PATH}/cmake/knime/binaries_mac.ini ${PAYLOAD}/binaries.ini
echo "done"

# get the search engines and add to bin/SEARCHENGINES
SEARCH_ARGS=1 
if [ $# -eq "$SEARCH_ARGS" ]; then
  SEARCH_ENGINES=$1
  if [ -e ${SEARCH_ENGINES}/OMSSA/omssacl -a -e ${SEARCH_ENGINES}/XTandem/tandem ]; then
    cp ${SEARCH_ENGINES}/OMSSA/omssacl ${SEARCH_ENGINES}/OMSSA/mods.xml ${SEARCH_ENGINES}/OMSSA/usermods.xml ${PAYLOAD}/bin/
    cp ${SEARCH_ENGINES}/XTandem/tandem ${PAYLOAD}/bin/
  else
    echo "No search engines found!!"
  fi
fi 

cd ${PAYLOAD}
zip -r binaries_mac_64.zip ./*
cd ${KNIME_PACKAGE}

echo -n "Cleaning payload .. "
rm -rf ${PAYLOAD}/binaries.ini
rm -rf ${PAYLOAD}/bin/
rm -rf ${PAYLOAD}/lib/
rm -rf ${PAYLOAD}/share/
echo "done"