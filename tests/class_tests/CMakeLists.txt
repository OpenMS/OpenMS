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
project("OpenMS_class_tests")

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
# set new CMAKE_RUNTIME_OUTPUT_DIRECTORY for tests and remember old setting
# Warning: in this directory the OpenMS.dll should already be present
#          (e.g. by using the Dll's POST_BUILD command to copy it)
#          otherwise the tests might not run!
set(_TMP_CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

#------------------------------------------------------------------------------
# Configure test file to get the TEST_DATA_PATH into the tests
set(CF_OPENMS_TEST_DATA_PATH "${PROJECT_SOURCE_DIR}/data/" CACHE INTERNAL "Path to the test data")
set (CONFIGURED_TEST_CONFIG_H "${PROJECT_BINARY_DIR}/include/OpenMS/test_config.h")
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/test_config.h.in ${CONFIGURED_TEST_CONFIG_H})

#------------------------------------------------------------------------------
# Copy dlls to test executables dir
if(MSVC)
  get_target_property(WIN32_DLLLOCATION OpenMS LOCATION)
  get_filename_component(WIN32_DLLPATH ${WIN32_DLLLOCATION} PATH)

  ## copy OpenMS.dll to test executables dir "$(TargetFileName)" is a placeholder filled by VS at runtime
  file(TO_NATIVE_PATH "${WIN32_DLLPATH}/$(TargetFileName)" DLL_SOURCE)
  file(TO_NATIVE_PATH "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$(ConfigurationName)/$(TargetFileName)" DLL_TARGET)

	## copy dlls to test runtime folder
	add_custom_command(TARGET OpenMS
                    POST_BUILD
										COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$(ConfigurationName)"
										COMMAND ${CMAKE_COMMAND} -E copy ${DLL_SOURCE} ${DLL_TARGET})
	add_custom_command(TARGET OpenMS_GUI
                    POST_BUILD
										COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$(ConfigurationName)"
										COMMAND ${CMAKE_COMMAND} -E copy ${DLL_SOURCE} ${DLL_TARGET})
	add_custom_command(TARGET OpenSwathAlgo
                    POST_BUILD
										COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$(ConfigurationName)"
                    COMMAND ${CMAKE_COMMAND} -E copy ${DLL_SOURCE} ${DLL_TARGET})
endif(MSVC)

#------------------------------------------------------------------------------
# get the test executables
include(executables.cmake)

#------------------------------------------------------------------------------
# required to build GUI tests (e.g., TOPPView_test)
add_definitions(-DQT_GUI_LIB)

#------------------------------------------------------------------------------
# Include directories for tests
set(OPENMS_CLASS_TESTS_INCLUDE_DIRECTORIES "${PROJECT_BINARY_DIR}/include/")
# add OpenMS directories
set(OPENMS_CLASS_TESTS_INCLUDE_DIRECTORIES "${OPENMS_CLASS_TESTS_INCLUDE_DIRECTORIES};${OPENMS_INCLUDE_DIRECTORIES}")
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
# VISUAL
add_custom_target(VISUAL_TEST)
add_dependencies(VISUAL_TEST ${visual_executables_list})

#------------------------------------------------------------------------------
# Add the actual tests
foreach(i ${TEST_executables})
	add_executable(${i} source/${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
	add_test(${i} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i})
	# only add OPENMP flags to gcc linker (execpt Mac OS X, due to compiler bug
	# see https://sourceforge.net/apps/trac/open-ms/ticket/280 for details)
  if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()
endforeach(i)

#------------------------------------------------------------------------------
# tests of classes from GUI lib
foreach(i ${visual_executables_list})
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
  target_include_directories(${i} PUBLIC ${OPENMS_GUI_INCLUDE_DIRECTORIES})
endforeach(i)

#------------------------------------------------------------------------------
# GUI tests special treatment - Apply MOC compiler to GUI
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Use QtTest lib
find_package(Qt4 REQUIRED QtTest)
if (NOT QT4_FOUND)
  message(STATUS "QtTest module not found!")
	message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -D QT_QMAKE_EXECUTABLE='<path_to_qmake(.exe)' <src-dir>")
endif()
include(${QT_USE_FILE})

#------------------------------------------------------------------------------
# Apply MOC compiler to GUI
set(GUI_executables_list_C)
foreach(i ${GUI_executables_list})
  list(APPEND GUI_executables_list_C GUI/${i}.C)
endforeach(i)

#------------------------------------------------------------------------------
## create mocced files - they are included manually in the implementation files
QT4_WRAP_CPP(mocced_GUI_executables_list_C ${GUI_executables_list_C} OPTIONS ${BOOST_MOC_ARGS})

foreach(i ${GUI_executables_list})
	SET_SOURCE_FILES_PROPERTIES(GUI/${i}.C PROPERTIES OBJECT_DEPENDS ${mocced_GUI_executables_list_C})
	add_executable(${i} source/GUI/${i}.C)
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES} ${QT_LIBRARIES})
  target_include_directories(${i} PUBLIC ${OPENMS_GUI_INCLUDE_DIRECTORIES})  
	add_test(${i} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i})
	# only add OPENMP flags to gcc linker (execpt Mac OS X, due to compiler bug
	# see https://sourceforge.net/apps/trac/open-ms/ticket/280 for details)
  if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()
endforeach(i)

#------------------------------------------------------------------------------
# GUI tests
add_custom_target(GUI_TEST)
add_dependencies(GUI_TEST ${GUI_executables_list})

#------------------------------------------------------------------------------
# add gui tests to executables list
set(TEST_executables ${TEST_executables} TOPPView_test)

#------------------------------------------------------------------------------
# restore old compiler flags
if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
	set(${CMAKE_CXX_FLAGS_RELEASE} ${_TMP_CMAKE_CXX_FLAGS_RELEASE})
endif()

#------------------------------------------------------------------------------
# "test_build" target
add_custom_target(test_build ALL)
add_dependencies(test_build ${TEST_executables} TOPP UTILS)

#------------------------------------------------------------------------------
# all target that we require to show up in the targets list
# and in the nightly build log (to test if they can be build properly)
add_custom_target(extra_targets ALL)
add_dependencies(extra_targets Tutorials_build)

#------------------------------------------------------------------------------
# add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${TEST_executables})
	list(APPEND sources_VS "${i}.C")
endforeach(i)
source_group("" FILES ${sources_VS})

#------------------------------------------------------------------------------
# add OpenSWATHAlgo specific tests
add_subdirectory(source/OPENSWATH)

#------------------------------------------------------------------------------
# restore old CMAKE_RUNTIME_OUTPUT_DIRECTORY
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${_TMP_CMAKE_RUNTIME_OUTPUT_DIRECTORY})
