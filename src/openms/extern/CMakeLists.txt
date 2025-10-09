# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Chris Bielow $
# $Authors: Stephan Aiche, Chris Bielow $
# --------------------------------------------------------------------------


project("OpenMS-3rd-party-libs")
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

# For use in the subdirs to cat their license
function(cat IN_FILE OUT_FILE)
    file(READ ${IN_FILE} CONTENTS)
    file(APPEND ${OUT_FILE} "${CONTENTS}")
endfunction()

# Prepare a temporary file to "cat" to:
set(TPLICENSEFILE ${OPENMS_HOST_BINARY_DIRECTORY}/thirdparty_licenses.txt)
file(WRITE ${TPLICENSEFILE} "The OpenMS library links statically to the following libraries with the following licenses:\n")
## TODO install it somewhere? Append to main license to show in installer?
##  At least we can use it in conda-forge after configuration

##
## internal packages (not checked if they exist externally)
##
add_subdirectory(evergreen/src) # evergreen Bayesian inference library (0.9)
add_subdirectory(eol-bspline)
add_subdirectory(IsoSpec)
add_subdirectory(GTE)
add_subdirectory(Quadtree)

if(ENABLE_TDL)
    add_subdirectory(tool_description_lib)
endif()

##
## external packages (with fallback option to local version)
##
option(USE_EXTERNAL_JSON "Use an external nlohmann::json library" OFF)
message(STATUS "Option USE_EXTERNAL_JSON: ${USE_EXTERNAL_JSON}")
if(USE_EXTERNAL_JSON)
    find_package(nlohmann_json REQUIRED GLOBAL)  ## global propagates the target up to the parent project
    message(STATUS "External 'NLohmann::json' at: '${nlohmann_json_DIR}'")
    target_compile_definitions(nlohmann_json::nlohmann_json INTERFACE JSON_USE_IMPLICIT_CONVERSIONS=0)
else()
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    add_subdirectory(nlohmann_json)
    set(JSON_USE_IMPLICIT_CONVERSIONS=0)
endif()

option(USE_EXTERNAL_SQLITECPP "Use an external SQLiteCpp library" OFF)
message(STATUS "Option USE_EXTERNAL_SQLITECPP: ${USE_EXTERNAL_SQLITECPP}")
if(USE_EXTERNAL_SQLITECPP)
    find_package(SQLiteCpp REQUIRED GLOBAL)      ## global propagates the target up to the parent project
    message(STATUS "External 'SQLiteCpp' at: '${SQLiteCpp_DIR}'")
else()
    add_subdirectory(SQLiteCpp)
endif()

option(USE_EXTERNAL_SIMDE "Use an external SIMDe library" OFF)
message(STATUS "Option USE_EXTERNAL_SIMDE: ${USE_EXTERNAL_SIMDE}")
if(USE_EXTERNAL_SIMDE)
    find_package(SIMDe REQUIRED GLOBAL)  ## global propagates the target up to the parent project
    message(STATUS "External 'SIMDe' at: '${SIMDe_INCLUDE_DIR}'")
else()
    add_subdirectory(simde)
endif()

option(ENABLE_ONNX "Enable support for the onnx runtime" OFF)
message(STATUS "Option ENABLE_ONNX: ${ENABLE_ONNX}")
if(ENABLE_ONNX)
    add_subdirectory(onnxruntime)
endif()
