#----------------------------------------------------------------------------------------------------------------------
# MIT License
#
# Copyright (c) 2021 Mark Schofield
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#----------------------------------------------------------------------------------------------------------------------
include_guard()

include("${CMAKE_CURRENT_LIST_DIR}/ToolchainCommon.cmake")

#[[====================================================================================================================
    toolchain_ensure_vswhere
    ------------------------

    Ensures that 'vswhere' is available and the VSWHERE_PATH variables is set in the cache. The following steps will be
    taken:

        1. If `VSWHERE_PATH` is set, then that value will be used.
        2. If 'vswhere' is found through `find_program`, that path will be used.
        3. 'vswhere' will be downloaded - specified by `VSWHERE_VERSION` and `VSWHERE_HASH` - into the
           `TOOLCHAIN_TOOLS_PATH` and `VSWHERE_PATH` will be set to the download location. If `TOOLCHAIN_TOOLS_PATH` is
           not set, then it will default to `${CMAKE_BINARY_DIR}/__tools`.

    Note: Since `CMAKE_BINARY_DIR` is platform specific, the default download location will change by platform,
    resulting in the tool being downloaded once for each platform that is built. Setting
    `TOOLCHAIN_TOOLS_PATH` to a platform-independent path (e.g. relative to the root of the repository) will
    allow 'vswhere' to be downloaded once for all platforms.
====================================================================================================================]]#
function(toolchain_ensure_vswhere)
    find_program(VSWHERE_PATH
            NAMES vswhere vswhere.exe
            HINTS "$ENV{ProgramFiles\(x86\)}/Microsoft Visual Studio/Installer"
            )

    # If vswhere isn't found, download it.
    #
    # vswhere.exe will be downloaded to `TOOLCHAIN_TOOLS_PATH`. If `TOOLCHAIN_TOOLS_PATH` is not set then it will be downloaded to '${CMAKE_BINARY_DIR}/__tools'.
    if(VSWHERE_PATH STREQUAL "VSWHERE_PATH-NOTFOUND")
        if(NOT VSWHERE_VERSION)
            set(VSWHERE_VERSION "2.8.4")
            set(VSWHERE_HASH "SHA256=E50A14767C27477F634A4C19709D35C27A72F541FB2BA5C3A446C80998A86419")
        else()
            if(NOT VSWHERE_HASH)
                message(FATAL_ERROR "VSWHERE_VERSION is set to ${VSWHERE_VERSION}. VSWHERE_HASH must be set if VSWHERE_VERSION is set.")
            endif()
        endif()

        if(NOT TOOLCHAIN_TOOLS_PATH)
            set(TOOLCHAIN_TOOLS_PATH "${CMAKE_BINARY_DIR}/__tools")
        endif()

        set(VSWHERE_PATH "${TOOLCHAIN_TOOLS_PATH}/vswhere.exe" CACHE FILEPATH "The location of 'vswhere.exe'" FORCE)

        toolchain_download_file(
                URL "https://github.com/microsoft/vswhere/releases/download/${VSWHERE_VERSION}/vswhere.exe"
                PATH ${VSWHERE_PATH}
                EXPECTED_HASH ${VSWHERE_HASH}
        )
    endif()
endfunction()

#[[====================================================================================================================
    findVisualStudio
    ----------------

    Finds a Visual Studio instance, and sets CMake variables based on properties of the found instance.

        findVisualStudio(
            [VERSION <version range>]
            [PRERELEASE <ON|OFF>]
            [PRODUCTS <products>]
            [REQUIRES <vs component>...]
            PROPERTIES
                <<vswhere property> <cmake variable>>
            )
====================================================================================================================]]#
function(findVisualStudio)
    set(OPTIONS)
    set(ONE_VALUE_KEYWORDS VERSION PRERELEASE PRODUCTS)
    set(MULTI_VALUE_KEYWORDS REQUIRES PROPERTIES)

    cmake_parse_arguments(PARSE_ARGV 0 FIND_VS "${OPTIONS}" "${ONE_VALUE_KEYWORDS}" "${MULTI_VALUE_KEYWORDS}")

    toolchain_ensure_vswhere()

    set(VSWHERE_COMMAND ${VSWHERE_PATH} -latest)

    if(FIND_VS_PRERELEASE)
        list(APPEND VSWHERE_COMMAND -prerelease)
    endif()

    if(FIND_VS_PRODUCTS)
        list(APPEND VSWHERE_COMMAND -products ${FIND_VS_PRODUCTS})
    endif()

    if(FIND_VS_REQUIRES)
        list(APPEND VSWHERE_COMMAND -requires ${FIND_VS_REQUIRES})
    endif()

    if(FIND_VS_VERSION)
        list(APPEND VSWHERE_COMMAND -version "${FIND_VS_VERSION}")
    endif()

    message(VERBOSE "findVisualStudio: VSWHERE_COMMAND = ${VSWHERE_COMMAND}")

    execute_process(
            COMMAND ${VSWHERE_COMMAND}
            OUTPUT_VARIABLE VSWHERE_OUTPUT
    )

    message(VERBOSE "findVisualStudio: VSWHERE_OUTPUT = ${VSWHERE_OUTPUT}")

    # Matches `VSWHERE_PROPERTY` in the `VSWHERE_OUTPUT` text in the format written by vswhere.
    # The matched value is assigned to the variable `VARIABLE_NAME` in the parent scope.
    function(getVSWhereProperty VSWHERE_OUTPUT VSWHERE_PROPERTY VARIABLE_NAME)
        string(REGEX MATCH "${VSWHERE_PROPERTY}: [^\r\n]*" VSWHERE_VALUE "${VSWHERE_OUTPUT}")
        string(REPLACE "${VSWHERE_PROPERTY}: " "" VSWHERE_VALUE "${VSWHERE_VALUE}")
        set(${VARIABLE_NAME} "${VSWHERE_VALUE}" PARENT_SCOPE)
    endfunction()

    while(FIND_VS_PROPERTIES)
        list(POP_FRONT FIND_VS_PROPERTIES VSWHERE_PROPERTY)
        list(POP_FRONT FIND_VS_PROPERTIES VSWHERE_CMAKE_VARIABLE)
        getVSWhereProperty("${VSWHERE_OUTPUT}" ${VSWHERE_PROPERTY} VSWHERE_VALUE)
        set(${VSWHERE_CMAKE_VARIABLE} ${VSWHERE_VALUE} PARENT_SCOPE)
    endwhile()
endfunction()