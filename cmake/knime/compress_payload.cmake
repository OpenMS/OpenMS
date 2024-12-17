# Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Stephan Aiche, Julianus Pfeuffer $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# include helper functions 
include ( ${SCRIPT_DIR}common.cmake )

set(required_variables "ARCH;PLATFORM;PAYLOAD_FOLDER")
check_variables(required_variables)

set(zip_file "${PAYLOAD_FOLDER}/binaries_${PLATFORM}_${ARCH}.zip")
file(TO_NATIVE_PATH ${zip_file} native_zip)

file(GLOB payload_content "${PAYLOAD_FOLDER}/*")
execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar "cf" "${native_zip}" --format=zip ${payload_content}
    WORKING_DIRECTORY "${PAYLOAD_FOLDER}"
)
# As always, CMake prematurely created an unusable new feature. Let's wait 5 more years until it is usable: https://gitlab.kitware.com/cmake/cmake/-/issues/21653
#file(ARCHIVE_CREATE OUTPUT ${native_zip} PATHS ${payload_content} FORMAT zip)

file(REMOVE_RECURSE ${payload_content})