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
# (OPENMS_WASM_THREADS, preset wasm32-emscripten-mt-release), or @openms/openms-wasm64-mt with
# threads and 64-bit memory (OPENMS_WASM_MEMORY64, preset wasm64-emscripten-mt-release): same
# files, same API.
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

set(_openms_wasm_sources
  bindings/binding_utils.cpp
  bindings/bind_chemistry.cpp
  bindings/bind_features.cpp
  bindings/bind_format.cpp
  bindings/bind_identification.cpp
  bindings/bind_kernel.cpp
  bindings/bind_misc.cpp
  bindings/bind_processing.cpp
  bindings/bind_prose.cpp
  bindings/bind_prose_tables.cpp
  bindings/js_types.cpp
)
add_executable(openms_wasm ${_openms_wasm_sources})
target_link_libraries(openms_wasm PRIVATE OpenMS)
openms_add_compiler_flags(openms_wasm)
openms_emscripten_link_malloc(openms_wasm) # the allocator of OPENMS_WASM_MALLOC (cmake/emscripten_malloc.cmake)

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
  # The virtual file system: WasmFS, the file system implemented in C++ inside the module, not
  # the JavaScript one. In the threaded build, the JavaScript file system lives on the main
  # thread and every file read of another thread waits for it; ProSE reads the spectra on a
  # helper thread while the main thread's team builds the index, and so got them too late to
  # restrict the index (a search against the human proteome took 8.7 s instead of 5.8 s). With
  # WasmFS, threads read directly. Directories of the host (Node.js) are mounted with
  # mountHostDirectory() (bindings/bind_misc.cpp); files of a browser are written to it.
  "-sWASMFS=1"
  "-sFORCE_FILESYSTEM=1"
  # Host files (mountHostDirectory()): an open or a close that fails reports its error with the sign
  # the backend expects. Emscripten's took the errno of an open for a file descriptor of the Node.js
  # process (wasmfs_node.js). SHELL: like the --pre-js below, so that CMake cannot drop it as a duplicate.
  "SHELL:--js-library \"${CMAKE_CURRENT_SOURCE_DIR}/wasmfs_node.js\""
  # FS.writeFile(), FS.write() and FS.read() throw ENOMEM where the memory cannot hold their buffer (Emscripten's
  # wrote the data from address 0 then), and FS.write() copies its data in one piece (wasmfs_fs.js)
  "SHELL:--js-library \"${CMAKE_CURRENT_SOURCE_DIR}/wasmfs_fs.js\""
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
set_property(TARGET openms_wasm APPEND PROPERTY LINK_DEPENDS "${_share_stamp}" "${CMAKE_CURRENT_SOURCE_DIR}/wasmfs_node.js"
             "${CMAKE_CURRENT_SOURCE_DIR}/wasmfs_fs.js")

## Threads: threads.js chooses the number of threads (Module.numThreads, OMP_NUM_THREADS)
## and the number of Web Workers started with the module (openmsThreadPoolSize): none
## in Node.js, which starts them on demand, and none on browser pages without shared
## memory (not cross-origin isolated), where the module then runs single-threaded. A
## fixed pool size would make the module hang at startup on such pages.
if(OPENMS_WASM_THREADS)
  ## ONNX Runtime starts a pool of threads of its own (onnx.js; before threads.js, which counts them).
  ## SHELL: CMake would drop the second "--pre-js" as a duplicate option
  if(WITH_ONNX)
    target_link_options(openms_wasm PRIVATE "SHELL:--pre-js \"${CMAKE_CURRENT_SOURCE_DIR}/onnx.js\"")
    set_property(TARGET openms_wasm APPEND PROPERTY LINK_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/onnx.js")
  endif()
  target_link_options(openms_wasm PRIVATE
    "--pre-js" "${CMAKE_CURRENT_SOURCE_DIR}/threads.js"
    "-sPTHREAD_POOL_SIZE=openmsThreadPoolSize"
    # threads.js starts the Workers of the threads from a Blob (Module.mainScriptUrlOrBlob), which
    # the module only reads if named here (Emscripten's default list plus that one)
    "-sINCOMING_MODULE_JS_API=ENVIRONMENT,arguments,canvas,dynamicLibraries,elementPointerLock,instantiateWasm,locateFile,monitorRunDependencies,noExitRuntime,noInitialRun,onAbort,onExit,onRuntimeInitialized,postRun,preInit,preRun,print,printErr,setStatus,statusMessage,stderr,stdin,stdout,thisProgram,wasm,websocket,mainScriptUrlOrBlob")
  set_property(TARGET openms_wasm APPEND PROPERTY LINK_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/threads.js")
  ## 64-bit memory (OPENMS_WASM_MEMORY64): a package of its own, so that an application can carry both. Memory64 needs
  ## Node.js 24 (Chrome 133, Firefox 134; not Safari).
  if(OPENMS_WASM_MEMORY64)
    set(OPENMS_WASM_PACKAGE_NAME "@openms/openms-wasm64-mt")
    set(OPENMS_WASM_PACKAGE_DESCRIPTION "OpenMS, the C++ library for computational mass spectrometry, compiled to multi-threaded WebAssembly (OpenMP) with 64-bit memory (Memory64, up to 16 GB) for web applications and Node.js")
    set(OPENMS_WASM_NODE_ENGINE ">=24")
  else()
    set(OPENMS_WASM_PACKAGE_NAME "@openms/openms-wasm-mt")
    set(OPENMS_WASM_PACKAGE_DESCRIPTION "OpenMS, the C++ library for computational mass spectrometry, compiled to multi-threaded WebAssembly (OpenMP) for web applications and Node.js")
  endif()
else()
  set(OPENMS_WASM_PACKAGE_NAME "@openms/openms-wasm")
  set(OPENMS_WASM_PACKAGE_DESCRIPTION "OpenMS, the C++ library for computational mass spectrometry, compiled to WebAssembly for web applications and Node.js")
endif()

#------------------------------------------------------------------------------
# One JavaScript API for wasm32 and wasm64: no integer of pointer size crosses to JavaScript, and no integer
# argument arrives unchecked
#
# long and unsigned long (size_t, ptrdiff_t) are 64 bits in the wasm64 build (OPENMS_WASM_MEMORY64), where
# Embind passes them as BigInt; the guards of the bound functions take and return them as numbers, and take
# every integer argument as number and check it (bindings/binding_utils.h), where Embind would truncate it.
# This check compiles the bindings once more (unoptimized, objects not linked) against a copy of Embind's wire.h
# in which the conversions of long and unsigned long are deleted, the conversions of the other integers (int,
# unsigned, short, char) from JavaScript to C++ (but not those of results), and memory views of long and
# unsigned long data (a BigUint64Array in wasm64, a Uint32Array in wasm32). Every other place where they would
# cross is a compile error, in both builds: e.g. a size_t as value of a val, an unguarded constructor<UInt>(),
# val::as<int>(), typed_memory_view() of size_t data. tests/check_dts.mjs (test openms_wasm_typescript) compares
# the declarations of the two builds.
option(OPENMS_WASM_CHECK_INTEGERS "Check at build time that no long or unsigned long (size_t) crosses from the OpenMS WebAssembly bindings to JavaScript, where wasm64 would pass it as BigInt, and that no integer argument from JavaScript bypasses the checks of the guards" ON)
if(OPENMS_WASM_CHECK_INTEGERS)
  set(_wire_h "${EMSCRIPTEN_SYSROOT}/include/emscripten/wire.h")
  set(_native_types [=[
EMSCRIPTEN_DEFINE_NATIVE_BINDING_TYPE(char);
EMSCRIPTEN_DEFINE_NATIVE_BINDING_TYPE(signed char);
EMSCRIPTEN_DEFINE_NATIVE_BINDING_TYPE(unsigned char);
EMSCRIPTEN_DEFINE_NATIVE_BINDING_TYPE(signed short);
EMSCRIPTEN_DEFINE_NATIVE_BINDING_TYPE(unsigned short);
EMSCRIPTEN_DEFINE_NATIVE_BINDING_TYPE(signed int);
EMSCRIPTEN_DEFINE_NATIVE_BINDING_TYPE(unsigned int);
EMSCRIPTEN_DEFINE_NATIVE_BINDING_TYPE(signed long);
EMSCRIPTEN_DEFINE_NATIVE_BINDING_TYPE(unsigned long);
]=])
  set(_checked_types [=[
// OpenMS (src/openms_wasm/CMakeLists.txt): long and unsigned long may not cross to JavaScript, where a
// wasm64 module passes them as BigInt, and no integer may come from JavaScript unchecked (Embind truncates it).
// Bind functions with OMS_FN, OMS_METHOD, OMS_CONSTRUCTOR, which take them as double and check them.
#define OPENMS_WASM_DELETED_BINDING_TYPE(type)                                 \
template<>                                                                     \
struct BindingType<type> {                                                     \
    typedef type WireType;                                                     \
    static WireType toWireType(const type& v, rvp::default_tag) = delete;      \
    static type fromWireType(WireType v) = delete;                             \
}
#define OPENMS_WASM_RESULT_BINDING_TYPE(type)                                  \
template<>                                                                     \
struct BindingType<type> {                                                     \
    typedef type WireType;                                                     \
    constexpr static WireType toWireType(const type& v, rvp::default_tag) {    \
        return v;                                                              \
    }                                                                          \
    static type fromWireType(WireType v) = delete;                             \
}
OPENMS_WASM_RESULT_BINDING_TYPE(char);
OPENMS_WASM_RESULT_BINDING_TYPE(signed char);
OPENMS_WASM_RESULT_BINDING_TYPE(unsigned char);
OPENMS_WASM_RESULT_BINDING_TYPE(signed short);
OPENMS_WASM_RESULT_BINDING_TYPE(unsigned short);
OPENMS_WASM_RESULT_BINDING_TYPE(signed int);
OPENMS_WASM_RESULT_BINDING_TYPE(unsigned int);
OPENMS_WASM_DELETED_BINDING_TYPE(signed long);
OPENMS_WASM_DELETED_BINDING_TYPE(unsigned long);
]=])
  set(_memory_view "constexpr bool typeSupportsMemoryView() {\n")
  set(_checked_memory_view "constexpr bool typeSupportsMemoryView() {\n    if (std::is_same<T, long>::value || std::is_same<T, unsigned long>::value) return false; // OpenMS: not of size_t data\n")
  file(READ "${_wire_h}" _wire)
  string(FIND "${_wire}" "${_native_types}" _at)
  string(FIND "${_wire}" "${_memory_view}" _at_view)
  if(_at EQUAL -1 OR _at_view EQUAL -1)
    message(WARNING "Emscripten's ${_wire_h} does not define the binding types of the integers or typeSupportsMemoryView() as expected: "
                    "the check of OPENMS_WASM_CHECK_INTEGERS is skipped (update it in src/openms_wasm/CMakeLists.txt).")
  else()
    string(REPLACE "${_native_types}" "${_checked_types}" _wire "${_wire}")
    string(REPLACE "${_memory_view}" "${_checked_memory_view}" _wire "${_wire}")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/integer_check/wire.h.new" "${_wire}")
    configure_file("${CMAKE_CURRENT_BINARY_DIR}/integer_check/wire.h.new"
                   "${CMAKE_CURRENT_BINARY_DIR}/integer_check/emscripten/wire.h" COPYONLY) # (rewritten if changed)
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_wire_h}")
    add_library(openms_wasm_integer_check OBJECT EXCLUDE_FROM_ALL ${_openms_wasm_sources} ${OPENMS_WASM_EXTRA_SOURCES})
    target_link_libraries(openms_wasm_integer_check PRIVATE OpenMS)
    openms_add_compiler_flags(openms_wasm_integer_check)
    target_include_directories(openms_wasm_integer_check BEFORE PRIVATE
      "${CMAKE_CURRENT_BINARY_DIR}/integer_check" "${CMAKE_CURRENT_SOURCE_DIR}/bindings")
    target_compile_options(openms_wasm_integer_check PRIVATE -O0 -g0)
    add_dependencies(openms_wasm openms_wasm_integer_check)
  endif()
endif()

#------------------------------------------------------------------------------
# The npm package around the module

if(NOT OPENMS_WASM_NODE_ENGINE)
  set(OPENMS_WASM_NODE_ENGINE ">=18")
endif()
configure_file(package.json.in "${OPENMS_WASM_PACKAGE_DIR}/package.json" @ONLY)
add_custom_command(TARGET openms_wasm POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/README.md" "${OPENMS_WASM_PACKAGE_DIR}/README.md"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${OPENMS_HOST_DIRECTORY}/LICENSE" "${OPENMS_WASM_PACKAGE_DIR}/LICENSE"
  # the demo web application (not part of the npm package: see "files" in package.json)
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/examples" "${OPENMS_WASM_PACKAGE_DIR}/examples"
  VERBATIM)
## ONNX Runtime: the PeptDeep models next to the module (models/, not part of the npm package
## either), for applications to load when needed: 22 MB that most uses of the module do not need.
if(WITH_ONNX)
  add_custom_command(TARGET openms_wasm POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${OPENMS_PEPTDEEP_MODEL_DIR}" "${OPENMS_WASM_PACKAGE_DIR}/models"
    VERBATIM)
endif()

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

  ## The TypeScript declarations (tests/check_dts.mjs): no BigInt, and the same as those of the build for the other
  ## memory (wasm32 for wasm64 and the other way round). OPENMS_WASM_COMPARE_DTS: the openms.d.ts of that build; by
  ## default the one of the build directory of the matching preset (build/wasm64-... for build/wasm32-...), compared
  ## if it is there.
  get_filename_component(_build_name "${OPENMS_HOST_BINARY_DIRECTORY}" NAME)
  get_filename_component(_build_parent "${OPENMS_HOST_BINARY_DIRECTORY}" DIRECTORY)
  set(_other_dts "")
  if(_build_name MATCHES "^wasm64-(.*)$")
    set(_other_dts "${_build_parent}/wasm32-${CMAKE_MATCH_1}/openms-wasm/openms.d.ts")
  elseif(_build_name MATCHES "^wasm32-(.*)$")
    set(_other_dts "${_build_parent}/wasm64-${CMAKE_MATCH_1}/openms-wasm/openms.d.ts")
  endif()
  set(OPENMS_WASM_COMPARE_DTS "" CACHE FILEPATH "openms.d.ts of the OpenMS WebAssembly build for the other memory (wasm32/wasm64), compared by the test openms_wasm_typescript (empty: the build directory of the matching preset, if built)")
  if(OPENMS_WASM_COMPARE_DTS)
    set(_compare "${OPENMS_WASM_COMPARE_DTS}")
  elseif(_other_dts)
    set(_compare "${_other_dts}" "--optional")
  else()
    set(_compare "")
  endif()
  add_test(NAME openms_wasm_typescript
    COMMAND "${OPENMS_WASM_NODE_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/check_dts.mjs" "${OPENMS_WASM_PACKAGE_DIR}/openms.d.ts" ${_compare})
endif()
