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

cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)
project("OpenMS_class_tests_openms")

#------------------------------------------------------------------------------
# -- BEGIN OF WORKAROUND --
# this is a workaround for a CMake-Bug (http://www.cmake.org/Bug/view.php?id=12419)
# , which prevents the orginal svn_revision_update target to be build prior to OpenMS.lib on Win-VS
# Delete when bug is fixed.
set(SVN_REVISION_FILE "${OPENMS_HOST_BINARY_DIRECTORY}/openms/include/OpenMS/openms_svn_revision.h")
add_custom_target(svn_revision_update_forTests ALL
	    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${OPENMS_HOST_DIRECTORY}
                               -DSVN_REVISION_FILE=${SVN_REVISION_FILE}
                               -DOPENMS_HAS_SVNVERSION=${OPENMS_HAS_SVNVERSION}
                               -DSVNVERSION_EXECUTABLE=${SVNVERSION_EXECUTABLE}
	                             -P ${OPENMS_HOST_DIRECTORY}/cmake/svnversion.cmake)
add_dependencies(OpenMS svn_revision_update_forTests)

# -- END OF WORKAROUND --
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Configure test file to get the TEST_DATA_PATH into the tests
set(CF_OPENMS_TEST_DATA_PATH "${PROJECT_SOURCE_DIR}/data/" CACHE INTERNAL "Path to the test data")
set (CONFIGURED_TEST_CONFIG_H "${PROJECT_BINARY_DIR}/include/OpenMS/test_config.h")
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/test_config.h.in ${CONFIGURED_TEST_CONFIG_H})

#------------------------------------------------------------------------------
# Some tests have a link dependency to boost so we find our own boost here
set(Boost_USE_STATIC_LIBS ${BOOST_USE_STATIC})
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)
add_definitions(/DBOOST_ALL_NO_LIB) ## disable auto-linking of boost libs (boost tends to guess wrong lib names)
set(Boost_COMPILER "")

# help boost finding it's packages
set(Boost_ADDITIONAL_VERSIONS "1.47.0" "1.48.0" "1.49.0" "1.50.0" "1.51.0" "1.52.0" "1.53.0" "1.54.0")

# 1st attempt does not explicitly requires boost to enable second check (see below)
find_package(Boost 1.42.0 COMPONENTS math_c99)

set(BOOST_MOC_ARGS "")

# see: https://bugreports.qt-project.org/browse/QTBUG-22829
# Confirmed only on mac os x and leads to problems on win32 and lnx
# so we handle it for now only on mac os x
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(BOOST_MOC_ARGS "-DBOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION")
endif()

#------------------------------------------------------------------------------
# get the test executables
include(executables.cmake)

#------------------------------------------------------------------------------
# Include directories for tests
set(OPENMS_CLASS_TESTS_INCLUDE_DIRECTORIES "${PROJECT_BINARY_DIR}/include/")
# add OpenMS directories
set(OPENMS_CLASS_TESTS_INCLUDE_DIRECTORIES "${OPENMS_CLASS_TESTS_INCLUDE_DIRECTORIES};${OPENMS_INCLUDE_DIRECTORIES};${Boost_INCLUDE_DIRS}")
include_directories(${OPENMS_CLASS_TESTS_INCLUDE_DIRECTORIES})

#------------------------------------------------------------------------------
# disable optimization for tests for gcc like compilers
if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
	set(_TMP_CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
	set(CMAKE_CXX_FLAGS_RELEASE "-O0")
endif()

#------------------------------------------------------------------------------
# CONCEPT
add_custom_target(CONCEPT_TEST)
add_dependencies(CONCEPT_TEST ${concept_executables_list})

#------------------------------------------------------------------------------
# DATASTRUCTURES
add_custom_target(DATASTRUCTURES_TEST)
add_dependencies(DATASTRUCTURES_TEST ${datastructures_executables_list})

#------------------------------------------------------------------------------
# METADATA
add_custom_target(METADATA_TEST)
add_dependencies(METADATA_TEST ${metadata_executables_list})

#------------------------------------------------------------------------------
# SYSTEM
add_custom_target(SYSTEM_TEST)
add_dependencies(SYSTEM_TEST ${system_executables_list})

#------------------------------------------------------------------------------
# KERNEL
add_custom_target(KERNEL_TEST)
add_dependencies(KERNEL_TEST ${kernel_executables_list})

#------------------------------------------------------------------------------
# FORMAT
add_custom_target(FORMAT_TEST)
add_dependencies(FORMAT_TEST ${format_executables_list})

#------------------------------------------------------------------------------
# MATH
add_custom_target(MATH_TEST)
add_dependencies(MATH_TEST ${math_executables_list})

#------------------------------------------------------------------------------
# FILTERING
add_custom_target(FILTERING_TEST)
add_dependencies(FILTERING_TEST ${filtering_executables_list})

#------------------------------------------------------------------------------
# COMPARISON
add_custom_target(COMPARISON_TEST)
add_dependencies(COMPARISON_TEST ${comparison_executables_list})

#------------------------------------------------------------------------------
# CHEMISTRY
add_custom_target(CHEMISTRY_TEST)
add_dependencies(CHEMISTRY_TEST ${chemistry_executables_list})

#------------------------------------------------------------------------------
# ANALYSIS
add_custom_target(ANALYSIS_TEST)
add_dependencies(ANALYSIS_TEST ${analysis_executables_list})

#------------------------------------------------------------------------------
# APPLICATIONS
add_custom_target(APPLICATIONS_TEST)
add_dependencies(APPLICATIONS_TEST ${applications_executables_list})

#------------------------------------------------------------------------------
# TRANSFORMATIONS
add_custom_target(TRANSFORMATIONS_TEST)
add_dependencies(TRANSFORMATIONS_TEST ${transformations_executables_list})

#------------------------------------------------------------------------------
# SIMULATION
add_custom_target(SIMULATION_TEST)
add_dependencies(SIMULATION_TEST ${simulation_executables_list})

#------------------------------------------------------------------------------
# Add the actual tests
foreach(_class_test ${TEST_executables})
	add_executable(${_class_test} source/${_class_test}.cpp)
	target_link_libraries(${_class_test} ${OPENMS_LIBRARIES})
	add_test(${_class_test} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_class_test})
	# only add OPENMP flags to gcc linker (execpt Mac OS X, due to compiler bug
	# see https://sourceforge.net/apps/trac/open-ms/ticket/280 for details)
  if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		set_target_properties(${_class_test} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()
endforeach(_class_test)

#------------------------------------------------------------------------------
# Tests with extra dependencies
#------------------------------------------------------------------------------
if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost was not found!")
endif()

#------------------------------------------------------------------------------
# 2 - add link dependency
foreach(t ${Boost_dependent_tests})
	target_link_libraries(${t} ${Boost_LIBRARIES})
endforeach()

#------------------------------------------------------------------------------
# restore old compiler flags
if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
	set(${CMAKE_CXX_FLAGS_RELEASE} ${_TMP_CMAKE_CXX_FLAGS_RELEASE})
endif()

#------------------------------------------------------------------------------
# all target that we require to show up in the targets list
# and in the nightly build log (to test if they can be build properly)
add_custom_target(extra_targets ALL)
add_dependencies(extra_targets Tutorials_build)

#------------------------------------------------------------------------------
# add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${TEST_executables})
	list(APPEND sources_VS "${i}.cpp")
endforeach(i)
source_group("" FILES ${sources_VS})
