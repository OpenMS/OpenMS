# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2024.
#
# This software is released under a three-clause BSD license.
# --------------------------------------------------------------------------
# pyOpenMS2: nanobind-based Python bindings for OpenMS
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20)

# Project name and version
project(pyOpenMS2
    VERSION 3.4.0
    LANGUAGES CXX
    DESCRIPTION "Python bindings for OpenMS using nanobind"
)

# C++17 is required for nanobind
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --------------------------------------------------------------------------
# Options
# --------------------------------------------------------------------------
option(PYOPENMS2_BUILD_TESTS "Build pyOpenMS2 tests" ON)
option(PYOPENMS2_GENERATE_STUBS "Generate .pyi stub files" ON)
set(PY_NUM_MODULES "8" CACHE STRING "Number of modules to split bindings into")

# --------------------------------------------------------------------------
# Find dependencies
# --------------------------------------------------------------------------
find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# Find or fetch nanobind
if(NOT TARGET nanobind::nanobind)
    include(FetchContent)
    FetchContent_Declare(
        nanobind
        GIT_REPOSITORY https://github.com/wjakob/nanobind.git
        GIT_TAG        v2.4.0
    )
    FetchContent_MakeAvailable(nanobind)
endif()

# Find OpenMS
find_package(OpenMS REQUIRED)

# Handle both namespaced and non-namespaced target names
if(TARGET OpenMS::OpenMS)
    set(OPENMS_TARGET OpenMS::OpenMS)
    set(OPENSWATHALGO_TARGET OpenMS::OpenSwathAlgo)
elseif(TARGET OpenMS)
    set(OPENMS_TARGET OpenMS)
    set(OPENSWATHALGO_TARGET OpenSwathAlgo)
else()
    message(FATAL_ERROR "OpenMS target not found")
endif()

# Note: Include directories are propagated automatically via target_link_libraries
# from the OpenMS targets' INTERFACE_INCLUDE_DIRECTORIES property

# --------------------------------------------------------------------------
# Generator: Parse .pxd files and generate nanobind C++ code
# --------------------------------------------------------------------------
set(PYOPENMS2_GENERATOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/generator")
set(PYOPENMS2_BINDINGS_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(PYOPENMS_PXD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../pyOpenMS/pxds")
set(PYOPENMS_ADDONS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../pyOpenMS/addons")

# Output sentinel file to track generation status
set(PYGEN_SENTINEL "${PYOPENMS2_BINDINGS_DIR}/.generated")

# Generated source files (one per module)
set(PYOPENMS2_GENERATED_SOURCES)
foreach(i RANGE 1 ${PY_NUM_MODULES})
    list(APPEND PYOPENMS2_GENERATED_SOURCES "${PYOPENMS2_BINDINGS_DIR}/module_${i}.cpp")
endforeach()
list(APPEND PYOPENMS2_GENERATED_SOURCES "${PYOPENMS2_BINDINGS_DIR}/main_module.cpp")

# Get OpenMS include directories for libclang parsing
get_target_property(OPENMS_INC_DIRS ${OPENMS_TARGET} INTERFACE_INCLUDE_DIRECTORIES)

# Build generator arguments (libclang is always used)
set(GENERATOR_ARGS
    --pxd-dir "${PYOPENMS_PXD_DIR}"
    --addons-dir "${PYOPENMS_ADDONS_DIR}"
    --output-dir "${PYOPENMS2_BINDINGS_DIR}"
    --num-modules ${PY_NUM_MODULES}
    --openms-include-dir ${OPENMS_INC_DIRS}
    --libclang-cache-dir "${CMAKE_CURRENT_BINARY_DIR}/libclang_cache"
    --all-classes
)
message(STATUS "pyOpenMS2: Using libclang for accurate C++ type parsing (with caching)")

# Collect .pxd and .pyx input files for dependency tracking
file(GLOB PYOPENMS_PXD_FILES "${PYOPENMS_PXD_DIR}/*.pxd")
file(GLOB PYOPENMS_ADDON_FILES "${PYOPENMS_ADDONS_DIR}/*.pyx")

# Custom command to run the generator
add_custom_command(
    OUTPUT ${PYGEN_SENTINEL} ${PYOPENMS2_GENERATED_SOURCES}
    COMMAND ${Python_EXECUTABLE} -m generator ${GENERATOR_ARGS}
    COMMAND ${CMAKE_COMMAND} -E touch ${PYGEN_SENTINEL}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS
        # Generator implementation files
        ${PYOPENMS2_GENERATOR_DIR}/pxd_parser.py
        ${PYOPENMS2_GENERATOR_DIR}/nanobind_emitter.py
        ${PYOPENMS2_GENERATOR_DIR}/cpp_parser.py
        ${PYOPENMS2_GENERATOR_DIR}/type_registry.py
        ${PYOPENMS2_GENERATOR_DIR}/addon_processor.py
        ${PYOPENMS2_GENERATOR_DIR}/__main__.py
        # Input files - regenerate when these change
        ${PYOPENMS_PXD_FILES}
        ${PYOPENMS_ADDON_FILES}
    COMMENT "Generating nanobind bindings from .pxd files"
    VERBATIM
)

# Custom target for generation
add_custom_target(pyopenms2_generate
    DEPENDS ${PYGEN_SENTINEL}
)

# --------------------------------------------------------------------------
# Type caster headers
# --------------------------------------------------------------------------
set(PYOPENMS2_TYPE_CASTERS
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/type_casters/string.h
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/type_casters/dposition.h
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/type_casters/datavalue.h
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/type_casters/stl_containers.h
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/type_casters/all_casters.h
)

# --------------------------------------------------------------------------
# Build split modules for parallel compilation and memory management
# --------------------------------------------------------------------------

# Each submodule is a separate shared library
foreach(i RANGE 1 ${PY_NUM_MODULES})
    set(MODULE_NAME "_pyopenms2_${i}")

    nanobind_add_module(${MODULE_NAME}
        NB_DOMAIN pyopenms  # Share types across modules
        ${PYOPENMS2_BINDINGS_DIR}/module_${i}.cpp
    )

    target_include_directories(${MODULE_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/bindings/type_casters
        ${PYOPENMS2_BINDINGS_DIR}
    )

    target_link_libraries(${MODULE_NAME} PRIVATE
        ${OPENMS_TARGET}
        ${OPENSWATHALGO_TARGET}
    )

    # Ensure generator runs first
    add_dependencies(${MODULE_NAME} pyopenms2_generate)

    # Set output directory
    set_target_properties(${MODULE_NAME} PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pyopenms"
    )

    # Install
    install(TARGETS ${MODULE_NAME}
        LIBRARY DESTINATION pyopenms
    )
endforeach()

# Main entry module that imports all submodules
nanobind_add_module(_pyopenms2
    NB_DOMAIN pyopenms
    ${PYOPENMS2_BINDINGS_DIR}/main_module.cpp
)

target_include_directories(_pyopenms2 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/type_casters
    ${PYOPENMS2_BINDINGS_DIR}
)

target_link_libraries(_pyopenms2 PRIVATE
    ${OPENMS_TARGET}
    ${OPENSWATHALGO_TARGET}
)

add_dependencies(_pyopenms2 pyopenms2_generate)

set_target_properties(_pyopenms2 PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pyopenms"
)

install(TARGETS _pyopenms2
    LIBRARY DESTINATION pyopenms
)

# --------------------------------------------------------------------------
# Python package files
# --------------------------------------------------------------------------
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/pyopenms/__init__.py
    ${CMAKE_CURRENT_BINARY_DIR}/pyopenms/__init__.py
    COPYONLY
)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/pyopenms/_dataframes.py
    ${CMAKE_CURRENT_BINARY_DIR}/pyopenms/_dataframes.py
    COPYONLY
)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/pyopenms/_dataframes_compat.py
    ${CMAKE_CURRENT_BINARY_DIR}/pyopenms/_dataframes_compat.py
    COPYONLY
)

# Copy addons
file(GLOB ADDON_FILES "${CMAKE_CURRENT_SOURCE_DIR}/pyopenms/addons/*.py")
foreach(addon_file ${ADDON_FILES})
    get_filename_component(addon_name ${addon_file} NAME)
    configure_file(
        ${addon_file}
        ${CMAKE_CURRENT_BINARY_DIR}/pyopenms/addons/${addon_name}
        COPYONLY
    )
endforeach()

# Install Python package files (including stubs and PEP 561 marker)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pyopenms/
    DESTINATION pyopenms
    FILES_MATCHING
        PATTERN "*.py"
        PATTERN "*.pyi"
        PATTERN "py.typed"
)

# --------------------------------------------------------------------------
# Stub generation
# --------------------------------------------------------------------------
if(PYOPENMS2_GENERATE_STUBS)
    # Collect all module targets as dependencies
    set(_STUB_DEPS _pyopenms2)
    foreach(i RANGE 1 ${PY_NUM_MODULES})
        list(APPEND _STUB_DEPS "_pyopenms2_${i}")
    endforeach()

    add_custom_target(pyopenms2_stubs
        COMMAND ${Python_EXECUTABLE} -m nanobind.stubgen
                -m pyopenms
                -r                          # recurse into _pyopenms2_N submodules
                -P                          # include private _pyopenms2_* members
                -O ${CMAKE_CURRENT_BINARY_DIR}/pyopenms  # stubs go alongside .so and .py files
                -M ${CMAKE_CURRENT_BINARY_DIR}/pyopenms/py.typed  # PEP 561 marker
                -p ${CMAKE_CURRENT_SOURCE_DIR}/stubs.patch       # fix singleton types etc.
        COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fix_stubs.py
                ${CMAKE_CURRENT_BINARY_DIR}/pyopenms   # post-process: fix keywords, NDArray types
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS ${_STUB_DEPS}
        COMMENT "Generating .pyi stub files for pyopenms"
    )
endif()

# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------
if(PYOPENMS2_BUILD_TESTS)
    enable_testing()

    # Find pytest
    execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import pytest"
        RESULT_VARIABLE PYTEST_NOT_FOUND
        OUTPUT_QUIET ERROR_QUIET
    )

    if(PYTEST_NOT_FOUND)
        message(WARNING "pytest not found, tests will not be available")
    else()
        add_test(
            NAME pyopenms2_pytest
            COMMAND ${CMAKE_COMMAND} -E env
                    PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}
                    ${Python_EXECUTABLE} -m pytest
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests -v
        )
    endif()
endif()

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
message(STATUS "")
message(STATUS "pyOpenMS2 Configuration Summary:")
message(STATUS "  Python:            ${Python_EXECUTABLE} (${Python_VERSION})")
message(STATUS "  NumPy:             ${Python_NumPy_INCLUDE_DIRS}")
message(STATUS "  nanobind:          ${nanobind_SOURCE_DIR}")
message(STATUS "  OpenMS:            ${OpenMS_DIR}")
message(STATUS "  Number of modules: ${PY_NUM_MODULES}")
message(STATUS "  Use Doxygen:       ${PYOPENMS2_USE_DOXYGEN}")
if(PYOPENMS2_USE_DOXYGEN)
    message(STATUS "  Doxygen XML dir:   ${PYOPENMS2_DOXYGEN_XML_DIR}")
endif()
message(STATUS "  Use libclang:      ${PYOPENMS2_USE_LIBCLANG}")
message(STATUS "  Build tests:       ${PYOPENMS2_BUILD_TESTS}")
message(STATUS "  Generate stubs:    ${PYOPENMS2_GENERATE_STUBS}")
message(STATUS "")
