# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

project("OpenMS_CLI")

#------------------------------------------------------------------------------
# The TOPP tool framework (TOPPBase, ToolHandler, INIUpdater, ParameterInformation,
# SearchEngineBase, TOPPExternalToolBase, MapAlignerBase, OpenSwathBase) as its own
# shared library. It links libOpenMS publicly and nothing in libOpenMS links back,
# so the scientific library's exported interface stays free of the command-line
# framework. TOPP tools, the GUI library and the documentation programs link
# OpenMS_CLI in-tree and OpenMS::OpenMS_CLI from an installed package.
#
# The include paths are unchanged: the headers stay under OpenMS/APPLICATIONS/.
# The export macro is OPENMS_CLI_DLLAPI from the generated
# OpenMS/APPLICATIONS/OpenMS_CLIConfig.h.

include(${PROJECT_SOURCE_DIR}/includes.cmake)

openms_add_library(TARGET_NAME OpenMS_CLI
                   SOURCE_FILES ${OpenMS_CLI_sources}
                   HEADER_FILES ${OpenMS_CLI_sources_h}
                   INTERNAL_INCLUDES ${PROJECT_SOURCE_DIR}/include
                                     ${PROJECT_BINARY_DIR}/include
                   LINK_LIBRARIES OpenMS
                   DLL_EXPORT_PATH "OpenMS/APPLICATIONS/")

# No QObject classes here; finding Qt for the GUI must not enable moc for this library.
set_target_properties(OpenMS_CLI PROPERTIES AUTOMOC OFF)

if (MSVC)
  ## same strictness as libOpenMS: unused parameters and locals are errors
  target_compile_options(OpenMS_CLI PRIVATE "/we4100" "/we4189")
endif()

openms_doc_path("${PROJECT_SOURCE_DIR}/include")
