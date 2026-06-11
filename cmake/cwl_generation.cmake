# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause

# Enables `make generate_cwl_files` to generate the files inside of `share/commonwl`
# and adds install targets for the cmake files
include(GNUInstallDirs)

# List all TOPP tools
set(executables ${TOPP_TOOLS})

# Try to find cwltool for validation (optional)
find_program(CWLTOOL_EXECUTABLE cwltool)
if(NOT CWLTOOL_EXECUTABLE)
  message(WARNING "cwltool not found. CWL file validation will be skipped.")
endif()

# Create a custom target
add_custom_target(generate_cwl_files DEPENDS TOPP)

# Create output folder (if not existing)
add_custom_command(TARGET generate_cwl_files POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory ${OPENMS_SHARE_DIR}/commonwl
)

# Walk through all tools and generate a CWL file
foreach(TOOL ${executables})
  # Add CWL generation
  add_custom_command(
    TARGET  generate_cwl_files POST_BUILD
    COMMAND ${OPENMS_BINARY_DIR}/${TOOL} -write_cwl ${OPENMS_SHARE_DIR}/commonwl
  )
  if(CWLTOOL_EXECUTABLE)
    add_custom_command(
      TARGET generate_cwl_files POST_BUILD
      # TODO: remove the no-warning flag
      COMMAND ${CWLTOOL_EXECUTABLE} --no-warning --validate ${OPENMS_SHARE_DIR}/commonwl/${TOOL}.cwl
      COMMENT "Validating ${OPENMS_SHARE_DIR}/commonwl/${TOOL}.cwl"
    )
  endif()

endforeach()
