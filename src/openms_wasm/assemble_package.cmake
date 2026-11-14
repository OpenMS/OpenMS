# Assembles the npm package @openms/openms-wasm from the packages of the two WebAssembly
# builds, single-threaded (preset wasm32-emscripten-release) and multi-threaded (preset
# wasm32-emscripten-mt-release), and the loader (loader/index.mjs), which starts the one
# the environment can run:
#
#   cmake -DSINGLE=build/wasm32-emscripten-release/openms-wasm
#         -DMULTI=build/wasm32-emscripten-mt-release/openms-wasm
#         -DOUTPUT=build/openms-wasm-package -P src/openms_wasm/assemble_package.cmake
#
# Layout: index.mjs, index.d.ts (loader), openms.d.ts (the module, the same for both builds),
# st/ and mt/ (openms.mjs, openms.wasm and, with OPENMS_WASM_SHARE=PRELOAD, openms.data),
# package.json, README.md, LICENSE.

cmake_minimum_required(VERSION 3.19)
foreach(_var SINGLE MULTI OUTPUT)
  if(NOT ${_var})
    message(FATAL_ERROR "assemble_package.cmake: set -D${_var}=<directory>")
  endif()
  get_filename_component(${_var} "${${_var}}" ABSOLUTE)
endforeach()

## the package.json of a build: its name and version
function(_openms_wasm_read_package dir out_name out_version)
  if(NOT EXISTS "${dir}/package.json" OR NOT EXISTS "${dir}/openms.mjs")
    message(FATAL_ERROR "${dir} is not the package directory of an OpenMS WebAssembly build (openms.mjs, package.json)")
  endif()
  file(READ "${dir}/package.json" _json)
  string(JSON _name GET "${_json}" name)
  string(JSON _version GET "${_json}" version)
  set(${out_name} "${_name}" PARENT_SCOPE)
  set(${out_version} "${_version}" PARENT_SCOPE)
endfunction()

_openms_wasm_read_package("${SINGLE}" _single_name _single_version)
_openms_wasm_read_package("${MULTI}" _multi_name _multi_version)
if(NOT _single_name STREQUAL "@openms/openms-wasm" OR NOT _multi_name STREQUAL "@openms/openms-wasm-mt")
  message(FATAL_ERROR "SINGLE has to be the single-threaded build (@openms/openms-wasm, is ${_single_name}), "
                      "MULTI the multi-threaded one (@openms/openms-wasm-mt, is ${_multi_name})")
endif()
if(NOT _single_version STREQUAL _multi_version)
  message(FATAL_ERROR "The builds are of different OpenMS versions: ${_single_version}, ${_multi_version}")
endif()

## OUTPUT is replaced: only if it is a previous output of this script (not, e.g., -DOUTPUT=.)
if(EXISTS "${OUTPUT}" AND NOT (EXISTS "${OUTPUT}/index.mjs" AND EXISTS "${OUTPUT}/package.json"))
  message(FATAL_ERROR "${OUTPUT} exists and is not an assembled package (index.mjs, package.json): refusing to delete it")
endif()
file(REMOVE_RECURSE "${OUTPUT}")
foreach(_build IN ITEMS st mt)
  if(_build STREQUAL "st")
    set(_dir "${SINGLE}")
  else()
    set(_dir "${MULTI}")
  endif()
  file(GLOB _files "${_dir}/openms.mjs" "${_dir}/openms.wasm" "${_dir}/openms.data")
  file(COPY ${_files} DESTINATION "${OUTPUT}/${_build}")
endforeach()
file(COPY "${SINGLE}/openms.d.ts" "${SINGLE}/README.md" "${SINGLE}/LICENSE"
          "${CMAKE_CURRENT_LIST_DIR}/loader/index.mjs" "${CMAKE_CURRENT_LIST_DIR}/loader/index.d.ts"
     DESTINATION "${OUTPUT}")
set(OPENMS_WASM_VERSION "${_single_version}")
configure_file("${CMAKE_CURRENT_LIST_DIR}/loader/package.json.in" "${OUTPUT}/package.json" @ONLY)
message(STATUS "Assembled ${_single_name} ${_single_version} in ${OUTPUT}")
