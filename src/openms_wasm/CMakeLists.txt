# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# OpenMS for JavaScript: the OpenMS library as WebAssembly module with Embind bindings,
# for web applications (and Node.js) that use OpenMS without native binaries.
# Built by the WebAssembly configuration (preset wasm32-emscripten-release,
# cmake/emscripten.cmake); see README.md.
#
# The target openms_wasm writes an npm package into <build>/openms-wasm:
#   openms.mjs    ES module, default export: the factory function createOpenMS()
#   openms.wasm   the compiled library
#                 (with OpenMS' shared data, share/OpenMS, at /share/OpenMS of its file system)
#   openms.d.ts   TypeScript declarations
#   package.json, README.md, LICENSE
# The package is @openms/openms-wasm, or @openms/openms-wasm-mt when built with threads
# (OPENMS_WASM_THREADS, preset wasm32-emscripten-mt-release): same files, same API.
#------------------------------------------------------------------------------

project("OpenMS_wasm")

set(OPENMS_WASM_PACKAGE_DIR "${OPENMS_HOST_BINARY_DIRECTORY}/openms-wasm"
    CACHE PATH "Output directory of the OpenMS npm package (WebAssembly module and JavaScript glue)")

#------------------------------------------------------------------------------
# Shared data: the parts of share/OpenMS the library reads at run time. Files only
# used by tools and the GUI are left out, and so are the HMDB mapping tables (49 MB),
# which only accurate mass search needs (OPENMS_WASM_SHARE_EXCLUDE removes them;
# AccurateMassSearchEngine users can mount their own files).
set(OPENMS_WASM_SHARE_DIRECTORIES CHEMISTRY CV MAPPING SCHEMAS
    CACHE STRING "Subdirectories of share/OpenMS preloaded into the module's file system")
set(OPENMS_WASM_SHARE_EXCLUDE "CHEMISTRY/HMDB2StructMapping.tsv;CHEMISTRY/HMDBMappingFile.tsv"
    CACHE STRING "Files of share/OpenMS (relative paths) not preloaded")
## How the module gets the data:
##  EMBED    inside openms.wasm: one self-contained file, which bundlers and CDNs handle
##           like any other WebAssembly file (default)
##  PRELOAD  a separate openms.data, downloaded when the module starts; it is looked up
##           relative to the page (or working directory): pass Module.locateFile
##  NONE     not packaged: write it to /share/OpenMS of the virtual file system yourself
##           before the library needs it
set(OPENMS_WASM_SHARE "EMBED" CACHE STRING "How the OpenMS WebAssembly module packages OpenMS' shared data: EMBED, PRELOAD or NONE")
set_property(CACHE OPENMS_WASM_SHARE PROPERTY STRINGS EMBED PRELOAD NONE)
if(NOT OPENMS_WASM_SHARE MATCHES "^(EMBED|PRELOAD|NONE)$")
  message(FATAL_ERROR "OPENMS_WASM_SHARE must be EMBED, PRELOAD or NONE (is '${OPENMS_WASM_SHARE}').")
endif()

set(_share_stage "${CMAKE_CURRENT_BINARY_DIR}/share/OpenMS")
set(_share_stamp "${CMAKE_CURRENT_BINARY_DIR}/share.stamp")
set(_share_inputs)
foreach(_dir IN LISTS OPENMS_WASM_SHARE_DIRECTORIES)
  file(GLOB_RECURSE _files CONFIGURE_DEPENDS "${OPENMS_HOST_DIRECTORY}/share/OpenMS/${_dir}/*")
  list(APPEND _share_inputs ${_files})
endforeach()
add_custom_command(OUTPUT "${_share_stamp}"
  COMMAND ${CMAKE_COMMAND}
    "-DSOURCE=${OPENMS_HOST_DIRECTORY}/share/OpenMS"
    "-DDESTINATION=${_share_stage}"
    "-DDIRECTORIES=${OPENMS_WASM_SHARE_DIRECTORIES}"
    "-DEXCLUDE=${OPENMS_WASM_SHARE_EXCLUDE}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/stage_share.cmake"
  COMMAND ${CMAKE_COMMAND} -E touch "${_share_stamp}"
  DEPENDS ${_share_inputs} "${CMAKE_CURRENT_SOURCE_DIR}/stage_share.cmake"
  COMMENT "Staging the shared data of the OpenMS WebAssembly module"
  VERBATIM)
add_custom_target(openms_wasm_share DEPENDS "${_share_stamp}")

#------------------------------------------------------------------------------
# The module

add_executable(openms_wasm
  bindings/binding_utils.cpp
  bindings/bind_chemistry.cpp
  bindings/bind_features.cpp
  bindings/bind_format.cpp
  bindings/bind_identification.cpp
  bindings/bind_kernel.cpp
  bindings/bind_misc.cpp
  bindings/bind_processing.cpp
  bindings/bind_prose.cpp
  bindings/js_types.cpp
)
target_link_libraries(openms_wasm PRIVATE OpenMS)
openms_add_compiler_flags(openms_wasm)

## Your own C++ code, compiled into the same module: sources with EMSCRIPTEN_BINDINGS
## blocks, which can use the helpers of bindings/binding_utils.h (see README.md).
set(OPENMS_WASM_EXTRA_SOURCES "" CACHE STRING "Additional C++ sources (absolute paths) compiled into the OpenMS WebAssembly module")
target_sources(openms_wasm PRIVATE ${OPENMS_WASM_EXTRA_SOURCES})
target_include_directories(openms_wasm PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/bindings")
add_dependencies(openms_wasm openms_wasm_share)

set_target_properties(openms_wasm PROPERTIES
  OUTPUT_NAME "openms"
  SUFFIX ".mjs"
  RUNTIME_OUTPUT_DIRECTORY "${OPENMS_WASM_PACKAGE_DIR}")

target_link_options(openms_wasm PRIVATE
  "-lembind"
  # an ES module exporting a factory: import createOpenMS from "openms.mjs"; const OpenMS = await createOpenMS();
  "-sMODULARIZE=1"
  "-sEXPORT_ES6=1"
  "-sEXPORT_NAME=createOpenMS"
  "-sENVIRONMENT=web,worker,node"
  # the virtual file system, with Node.js directories (NODEFS) and browser File objects (WORKERFS) as mountable backends
  "-sFORCE_FILESYSTEM=1"
  "-lnodefs.js"
  "-lworkerfs.js"
  "-sEXPORTED_RUNTIME_METHODS=FS,ENV"
  # A library module: keep the runtime alive between calls.
  "-sEXIT_RUNTIME=0"
  "--emit-tsd" "${OPENMS_WASM_PACKAGE_DIR}/openms.d.ts"
)
if(OPENMS_WASM_SHARE STREQUAL "EMBED")
  target_link_options(openms_wasm PRIVATE "--embed-file" "${_share_stage}@/share/OpenMS")
elseif(OPENMS_WASM_SHARE STREQUAL "PRELOAD")
  target_link_options(openms_wasm PRIVATE "--preload-file" "${_share_stage}@/share/OpenMS")
endif()
# Relink when the data or the link options change.
set_property(TARGET openms_wasm APPEND PROPERTY LINK_DEPENDS "${_share_stamp}")

## Threads: threads.js chooses the number of threads (Module.numThreads, OMP_NUM_THREADS)
## and the number of Web Workers started with the module (openmsThreadPoolSize): none
## in Node.js, which starts them on demand, and none on browser pages without shared
## memory (not cross-origin isolated), where the module then runs single-threaded. A
## fixed pool size would make the module hang at startup on such pages.
if(OPENMS_WASM_THREADS)
  target_link_options(openms_wasm PRIVATE
    "--pre-js" "${CMAKE_CURRENT_SOURCE_DIR}/threads.js"
    "-sPTHREAD_POOL_SIZE=openmsThreadPoolSize")
  set_property(TARGET openms_wasm APPEND PROPERTY LINK_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/threads.js")
  set(OPENMS_WASM_PACKAGE_NAME "@openms/openms-wasm-mt")
  set(OPENMS_WASM_PACKAGE_DESCRIPTION "OpenMS, the C++ library for computational mass spectrometry, compiled to multi-threaded WebAssembly (OpenMP) for web applications and Node.js")
else()
  set(OPENMS_WASM_PACKAGE_NAME "@openms/openms-wasm")
  set(OPENMS_WASM_PACKAGE_DESCRIPTION "OpenMS, the C++ library for computational mass spectrometry, compiled to WebAssembly for web applications and Node.js")
endif()

#------------------------------------------------------------------------------
# The npm package around the module

configure_file(package.json.in "${OPENMS_WASM_PACKAGE_DIR}/package.json" @ONLY)
add_custom_command(TARGET openms_wasm POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/README.md" "${OPENMS_WASM_PACKAGE_DIR}/README.md"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${OPENMS_HOST_DIRECTORY}/LICENSE" "${OPENMS_WASM_PACKAGE_DIR}/LICENSE"
  # the demo web application (not part of the npm package: see "files" in package.json)
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/examples" "${OPENMS_WASM_PACKAGE_DIR}/examples"
  VERBATIM)

#------------------------------------------------------------------------------
# Tests: the bindings, run with Node.js (node --test)

## HINTS are directories; EMSDK_NODE (emsdk_env.sh) and NODE_JS_EXECUTABLE (Emscripten's
## toolchain file) are paths of the node binary.
set(_openms_wasm_node_hints)
foreach(_node IN ITEMS "$ENV{EMSDK_NODE}" "${NODE_JS_EXECUTABLE}")
  if(_node)
    get_filename_component(_node_dir "${_node}" DIRECTORY)
    list(APPEND _openms_wasm_node_hints "${_node_dir}")
  endif()
endforeach()
find_program(OPENMS_WASM_NODE_EXECUTABLE NAMES node
  HINTS ${_openms_wasm_node_hints}
  DOC "Node.js for the tests of the OpenMS WebAssembly module")
if(OPENMS_WASM_NODE_EXECUTABLE)
  file(GLOB _openms_wasm_tests CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.test.mjs")
  add_test(NAME openms_wasm_bindings
    COMMAND "${OPENMS_WASM_NODE_EXECUTABLE}" --test --test-reporter=spec ${_openms_wasm_tests}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  set_tests_properties(openms_wasm_bindings PROPERTIES
    ENVIRONMENT "OPENMS_WASM_MODULE=${OPENMS_WASM_PACKAGE_DIR}/openms.mjs;OPENMS_TEST_DATA=${OPENMS_HOST_DIRECTORY}/src/tests/class_tests/openms/data"
    TIMEOUT 600)
endif()
