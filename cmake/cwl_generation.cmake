# Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause

# Enables `make generate_cwl_files` to generate the files inside of `workflow/cwl`

# List all TOPP tools
set(executables ${TOPP_TOOLS})

# Tools that can't export proper CWL files
list(REMOVE_ITEM executables GenericWrapper)


# Create a custom target
add_custom_target(generate_cwl_files DEPENDS TOPP)

# Create output folder (if not existing)
add_custom_command(TARGET generate_cwl_files POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/workflow/cwl
)

# Walk through all tools and generate a CWL file
foreach(TOOL ${executables})
  add_custom_command(
    TARGET  generate_cwl_files POST_BUILD
    COMMAND ${OPENMS_BINARY_DIR}/${TOOL} -write_cwl ${CMAKE_CURRENT_SOURCE_DIR}/workflow/cwl
  )
endforeach()

