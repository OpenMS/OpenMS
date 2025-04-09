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

#[[====================================================================================================================
    toolchain_update_file
    ---------------------
    Updates the file at the given path to have the given contents. If the contents already match, then the file is not
    touched.

        toolchain_update_file(
            PATH <file path>
            CONTENT <content>
        )
====================================================================================================================]]#
function(toolchain_update_file)
    set(OPTIONS)
    set(ONE_VALUE_KEYWORDS CONTENT PATH)
    set(MULTI_VALUE_KEYWORDS)

    cmake_parse_arguments(PARSE_ARGV 0 UPDATE_FILE "${OPTIONS}" "${ONE_VALUE_KEYWORDS}" "${MULTI_VALUE_KEYWORDS}")

    if(NOT (EXISTS ${UPDATE_FILE_PATH}))
        message(VERBOSE "update_file: Creating ${UPDATE_FILE_PATH}")
        file(WRITE ${UPDATE_FILE_PATH} ${UPDATE_FILE_CONTENT})
    else()
        file(READ ${UPDATE_FILE_PATH} CURRENT_CONTENT)
        if(NOT (CURRENT_CONTENT STREQUAL UPDATE_FILE_CONTENT))
            message(VERBOSE "update_file: Updating ${UPDATE_FILE_PATH}")
            file(WRITE ${UPDATE_FILE_PATH} ${UPDATE_FILE_CONTENT})
        endif()
    endif()
endfunction()

#[[====================================================================================================================
    toolchain_find_powershell
    -------------------------
    Searches for PowerShell.

        toolchain_find_powershell(
            <VARIABLE_NAME>
        )

    If the POWERSHELL_PATH value is set, that will be preferred. The VARIABLE_NAME will be set in the parent scope
    with the result.
====================================================================================================================]]#
function(toolchain_find_powershell VARIABLE_NAME)
    if(NOT EXISTS ${POWERSHELL_PATH})
        find_program(
                POWERSHELL_PATH
                NAMES pwsh.exe powershell.exe
                DOC "The path to PowerShell"
                REQUIRED
        )
    endif()
    message(VERBOSE "toolchain_find_powershell: Using PowerShell from: ${POWERSHELL_PATH}")
    set(${VARIABLE_NAME} ${POWERSHELL_PATH} PARENT_SCOPE)
endfunction()

#[[====================================================================================================================
    toolchain_download_file
    -----------------------
    Downloads a file to the given path.

        toolchain_download_file(
            URL <url>
            PATH <file path>
            EXPECTED_HASH <algorithm>=<hash>
        )

    The URL and EXPECTED_HASH values will be written to "<file path>.key" forcing a re-download if the values change.
====================================================================================================================]]#
function(toolchain_download_file)
    set(OPTIONS)
    set(ONE_VALUE_KEYWORDS URL PATH EXPECTED_HASH)
    set(MULTI_VALUE_KEYWORDS)

    cmake_parse_arguments(PARSE_ARGV 0 DOWNLOAD_FILE "${OPTIONS}" "${ONE_VALUE_KEYWORDS}" "${MULTI_VALUE_KEYWORDS}")

    toolchain_update_file(
            PATH "${DOWNLOAD_FILE_PATH}.key"
            CONTENT "URL ${DOWNLOAD_FILE_URL}:EXPECTED_HASH ${DOWNLOAD_FILE_EXPECTED_HASH}"
    )

    if("${DOWNLOAD_FILE_PATH}.key" IS_NEWER_THAN ${DOWNLOAD_FILE_PATH})
        message(VERBOSE "download_file: ${DOWNLOAD_FILE_URL} --> ${DOWNLOAD_FILE_PATH}")
        file(REMOVE ${DOWNLOAD_FILE_PATH})
        file(DOWNLOAD
                ${DOWNLOAD_FILE_URL}
                ${DOWNLOAD_FILE_PATH}
                EXPECTED_HASH ${DOWNLOAD_FILE_EXPECTED_HASH}
                )
    endif()
endfunction()