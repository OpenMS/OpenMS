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
# $Maintainer: Witold Wolski $
# $Authors: Witold Wolski, Stephan Aiche $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.8.0 FATAL_ERROR)
project("OpenMS_class_tests_superhirn")

#------------------------------------------------------------------------------
# Configure test file to get the TEST_DATA_PATH into the tests
set(CF_SUPERHIRN_TEST_DATA_PATH "${PROJECT_SOURCE_DIR}/data/")
set (CONFIGURED_TEST_CONFIG_H "${PROJECT_BINARY_DIR}/test_config.h")
configure_file(${PROJECT_SOURCE_DIR}/test_config.h.in ${CONFIGURED_TEST_CONFIG_H})

# --------------------------------------------------------------------------
# add include directories
include_directories(SYSTEM ${OpenMS_INCLUDE_DIRECTORIES}
                           ${SuperHirn_INCLUDE_DIRECTORIES}
                           ${PROJECT_BINARY_DIR})

# --------------------------------------------------------------------------
set(superhirn_tests
  # superhirn classes
  SuperhirnBackgroundControl_test
  SuperhirnBackgroundIntensityBin_test
  SuperhirnCentroidData_test
  SuperhirnCentroidPeak_test
  SuperhirnClusteredMS2ConsensusSpectrum_test
  SuperhirnDeisotoper_test
  SuperhirnFT_PeakDetectController_test
  SuperhirnIsotopicDist_test
  SuperhirnLCMSCData_test
  SuperhirnLC_MS_test
  SuperhirnLC_elution_peak_test
  SuperhirnMS1_feature_merger_test
  SuperhirnMS2ConsensusSpectrum_test
  SuperhirnMS2Fragment_test
  SuperhirnMS2_feature_test
  SuperhirnMS2Info_test
  SuperhirnProcess_Data_test
  SuperhirnRawData_test
  SuperhirnconsensIsotopePattern_test
  SuperhirnfeatureLCprofile_test
  Superhirnfeature_test
  Superhirnms_peak_test
  # superhirn featurefinder
  Superhirn_test
)

#------------------------------------------------------------------------------
# QT dependencies
find_package(Qt5 COMPONENTS Core Network REQUIRED)
if (NOT Qt5Network_FOUND)
  message(STATUS "QtNetwork module not found!")
  message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -D QT_QMAKE_EXECUTABLE='<path_to_qmake(.exe)' <src-dir>")
endif()

# --------------------------------------------------------------------------
# add targets for the executables
foreach(i ${superhirn_tests})
  add_executable(${i} source/${i}.cpp)
  target_link_libraries(${i} SuperHirn OpenMS)
  target_compile_definitions(${i} PRIVATE BOOST_ALL_NO_LIB)
  add_test(${i} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i})
endforeach(i)
