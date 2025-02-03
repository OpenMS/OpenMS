# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright OpenMS Inc. -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-present.
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

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_class_tests_openms")

#------------------------------------------------------------------------------
# Configure test file to get the TEST_DATA_PATH into the tests
set(CF_OPENMS_TEST_DATA_PATH "${PROJECT_SOURCE_DIR}/data/")
set (CONFIGURED_TEST_CONFIG_H "${PROJECT_BINARY_DIR}/include/OpenMS/test_config.h")
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/test_config.h.in ${CONFIGURED_TEST_CONFIG_H})

#------------------------------------------------------------------------------
# get the test executables
include(executables.cmake)

#------------------------------------------------------------------------------
# Include directories for tests
set(OPENMS_CLASS_TESTS_INTERNAL_INCLUDE_DIRECTORIES "${PROJECT_BINARY_DIR}/include/")
# add OpenMS directories
set(OPENMS_CLASS_TESTS_EXTERNAL_INCLUDE_DIRECTORIES "${OpenMS_INCLUDE_DIRECTORIES};${Boost_INCLUDE_DIRS}")
include_directories(${OPENMS_CLASS_TESTS_INTERNAL_INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${OPENMS_CLASS_TESTS_EXTERNAL_INCLUDE_DIRECTORIES})

#------------------------------------------------------------------------------
# disable optimization for tests for gcc like compilers
if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
  set(_TMP_CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
  set(CMAKE_CXX_FLAGS_RELEASE "-O0")
endif()

#------------------------------------------------------------------------------
# Add the actual tests
foreach(_class_test ${TEST_executables})
  add_executable(${_class_test} source/${_class_test}.cpp)
  target_link_libraries(${_class_test} ${OpenMS_LIBRARIES})
  add_test(${_class_test} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_class_test})
  # only add OPENMP flags to gcc linker (except Mac OS X, due to compiler bug
  # see https://sourceforge.net/apps/trac/open-ms/ticket/280 for details)
  if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${_class_test} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endforeach(_class_test)

#------------------------------------------------------------------------------
# some tests need special treatment

# StopWatch_test needs to be run alone to avoid timing errors
set_tests_properties(StopWatch_test PROPERTIES RUN_SERIAL 1)

#------------------------------------------------------------------------------
# 2 - add link dependencies (TODO could be more finegrained for boost)
foreach(t ${Boost_dependent_tests})
  foreach(comp ${OpenMS_BOOST_COMPONENTS})
    target_link_libraries(${t} Boost::${comp})
  endforeach()
endforeach()

target_link_libraries(Base64_test ZLIB::ZLIB)
target_link_libraries(LPWrapper_test CoinOR::CoinOR)
target_link_libraries(SpectraSTSimilarityScore_test Eigen3::Eigen)
target_link_libraries(BinnedSpectrum_test Eigen3::Eigen)

if (WITH_HDF5)
  target_link_libraries(HDF5_test HDF5::HDF5)
endif()

#------------------------------------------------------------------------------
# restore old compiler flags
if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
  set(${CMAKE_CXX_FLAGS_RELEASE} ${_TMP_CMAKE_CXX_FLAGS_RELEASE})
endif()

#------------------------------------------------------------------------------
# add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${TEST_executables})
  list(APPEND sources_VS "${i}.cpp")
endforeach(i)
source_group("" FILES ${sources_VS})
