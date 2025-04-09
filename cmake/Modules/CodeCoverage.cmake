# Copyright (c) 2012 - 2015, Lars Bilke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#
# 2012-01-31, Lars Bilke
# - Enable Code Coverage
#
# 2013-09-17, Joakim Söderberg
# - Added support for Clang.
# - Some additional usage instructions.
#
# 2016-10-11, Julianus Pfeuffer
# - Adaption and restructuring of dependencies to fit OpenMS
#
# USAGE:

# 0. (Mac only) If you use Xcode 5.1 make sure to patch geninfo as described here:
#      http://stackoverflow.com/a/22404544/80480
#
# 1. Copy this file into your cmake modules path.
#
# 2. Add the following line to your CMakeLists.txt:
#      INCLUDE(CodeCoverage)
#
# 3. Set compiler flags to turn off optimization and enable coverage
#   (also use Debug build type to ensure the first two arguments):
#
#    SET(CMAKE_CXX_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
#    SET(CMAKE_C_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
#
# 3. Use the function SETUP_TARGET_FOR_COVERAGE to create a custom make target
#    which runs your test executable and produces a lcov code coverage report:
#    Example:
#    SETUP_TARGET_FOR_COVERAGE(
#               my_coverage_target      # Name for custom target.
#               coverage                # Name of output directory.
#               /home/my_external_libs  # Semicolon seperated paths to exclude external sources from report.
#               )
#
# 4. Build a Debug build and run tests:
#    cmake -DCMAKE_BUILD_TYPE=Debug ..
#    make
#    make test
#    make my_coverage_target
#
#

# Check prereqs
FIND_PROGRAM( GCOV_PATH gcov )
FIND_PROGRAM( LCOV_PATH lcov )
FIND_PROGRAM( GENHTML_PATH genhtml )

IF(NOT GCOV_PATH)
    MESSAGE(FATAL_ERROR "gcov not found! Aborting...")
ENDIF() # NOT GCOV_PATH

IF(NOT CMAKE_COMPILER_IS_GNUCXX)
    IF(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS 3.0.0)
        IF(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR "${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS 5.1.0)
          MESSAGE(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} is neither GNU gcc nor Clang Version > 3.0.0. No support for gcov Coverage analysis.")
        ENDIF()
    ENDIF()
ENDIF() # NOT CMAKE_COMPILER_IS_GNUCXX

SET(CMAKE_CXX_FLAGS_COVERAGE
    "-g -O0 --coverage -fprofile-arcs -ftest-coverage"
    CACHE STRING "Flags used by the C++ compiler during coverage builds."
    FORCE )
SET(CMAKE_C_FLAGS_COVERAGE
    "-g -O0 --coverage -fprofile-arcs -ftest-coverage"
    CACHE STRING "Flags used by the C compiler during coverage builds."
    FORCE )
SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE
    ""
    CACHE STRING "Flags used for linking binaries during coverage builds."
    FORCE )
SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE
    ""
    CACHE STRING "Flags used by the shared libraries linker during coverage builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE )

IF ( NOT (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "Coverage"))
  MESSAGE( WARNING "Code coverage results with an optimized (non-Debug) build may be misleading" )
ENDIF() # NOT CMAKE_BUILD_TYPE STREQUAL "Debug"


# Param _targetname         The name of new the custom make target
# Param _outputname         lcov output is generated as _outputname.info
#                           HTML report is generated in _outputname/index.html
# Param _addignorelibpaths  Adds these paths to the --remove option of lcov to e.g. exclude external sources

FUNCTION(SETUP_TARGET_FOR_COVERAGE _targetname _outputname _addignorelibpaths)

    IF(NOT LCOV_PATH)
        MESSAGE(FATAL_ERROR "lcov not found! Aborting...")
    ENDIF() # NOT LCOV_PATH

    IF(NOT GENHTML_PATH)
        MESSAGE(FATAL_ERROR "genhtml not found! Aborting...")
    ENDIF() # NOT GENHTML_PATH

    SET(coverage_info "${CMAKE_BINARY_DIR}/${_outputname}.info")
    SET(coverage_cleaned "${coverage_info}.cleaned")

    SET(ignorelibpaths \"tests/*\" \"/usr/*\" \"/Applications/*\")
    foreach(libpath ${_addignorelibpaths})
        list(APPEND ignorelibpaths \"${libpath}/*\")
    endforeach()

    ## Workaround that CMake does not complain during configure, that the log is missing.
    SET_SOURCE_FILES_PROPERTIES(
        ${CMAKE_CURRENT_BINARY_DIR}/Testing/Temporary/LastTest.log
        PROPERTIES GENERATED TRUE
    )

    ADD_CUSTOM_COMMAND(
        OUTPUT ${_outputname}/index.html
        
        # Capturing lcov counters and generating report
        COMMAND ${LCOV_PATH} --directory . --capture --output-file ${coverage_info}
        # Removing external sources
        COMMAND ${LCOV_PATH} --remove ${coverage_info} ${ignorelibpaths} --output-file ${coverage_cleaned}
        # Generating html
        COMMAND ${GENHTML_PATH} -o ${_outputname} ${coverage_cleaned}
        # Remove temporaries
        COMMAND ${CMAKE_COMMAND} -E remove ${coverage_info} ${coverage_cleaned}
        # Compares the timestamps of the last tests with the index.html of the generated report and
        # only rebuilds if tests were performed after last report generation.
        MAIN_DEPENDENCY ${CMAKE_CURRENT_BINARY_DIR}/Testing/Temporary/LastTest.log
        VERBATIM
        COMMENT "Coverage data outdated. Processing code coverage counters and generating report."
    )

    # Setup target
    ADD_CUSTOM_TARGET(${_targetname}
        # Depends on the output of the previous command. Always checks if this custom_command needs to be re-executed.
        DEPENDS testsExecuted ${_outputname}/index.html

        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Coverage report up-to-date. Re-run tests if you need a new report."
    )
    
    # This target is basically there for a better error message.
    # Otherwise you get "No target for Testing/Temporary/LastTest.log"
    # Workaround because you can not depend on internal targets like "test" (https://gitlab.kitware.com/cmake/cmake/issues/8438)
    # Alternative: Always auto-invoke "make test" before. But in some scenarios you already have them built already.
    ADD_CUSTOM_TARGET(testsExecuted
        COMMAND ${CMAKE_COMMAND} -E md5sum "${CMAKE_CURRENT_BINARY_DIR}/Testing/Temporary/LastTest.log"
        COMMENT "Checking existence of test timestamp. If this step fails, please run 'make test' (again)."
    )

    # Show info where to find the report and clean up.
    ADD_CUSTOM_COMMAND(TARGET ${_targetname} POST_BUILD
        COMMAND ${LCOV_PATH} --directory . --zerocounters
        COMMENT "Cleaning up counters.. for another fresh coverage scan please execute make test again.\nOpen ./${_outputname}/index.html in your browser to view the coverage report."
    )

ENDFUNCTION() # SETUP_TARGET_FOR_COVERAGE
