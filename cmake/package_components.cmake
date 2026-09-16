# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------

cpack_add_install_type(recommended DISPLAY_NAME "Recommended")
cpack_add_install_type(full DISPLAY_NAME "Full")
cpack_add_install_type(minimal DISPLAY_NAME "Minimal")

## TODO group some components like "OpenMS_header", "OpenSWATH_header", "thirdparty_headers"...
## TODO do components more fine-grained, so you can install only OpenSwath parts, or only TOPPView but not TOPPAS, etc
cpack_add_component(share
                DISPLAY_NAME "OpenMS shared files"
                DESCRIPTION "OpenMS shared files"
                INSTALL_TYPES recommended full minimal
                )
cpack_add_component(library
                DISPLAY_NAME "Libraries"
                DESCRIPTION "The OpenMS core libraries"
                INSTALL_TYPES recommended full minimal
                )
## the layers above the core library (see cmake/install_macros.cmake)
cpack_add_component(library_cli
                DISPLAY_NAME "TOPP tool framework library"
                DESCRIPTION "The TOPP tool framework library (OpenMS_CLI), needed by the TOPP tools"
                DEPENDS library
                INSTALL_TYPES recommended full minimal
                )
if(WITH_GUI)
  cpack_add_component(library_gui
                  DISPLAY_NAME "GUI library"
                  DESCRIPTION "The GUI library (OpenMS_GUI), needed by TOPPView, TOPPAS and the other GUI applications"
                  DEPENDS library_cli
                  INSTALL_TYPES recommended full minimal
                  )
endif()
cpack_add_component(applications
                DISPLAY_NAME "OpenMS binaries"
                DESCRIPTION "OpenMS binaries including TOPP tools, TOPPView and TOPPAS."
                INSTALL_TYPES recommended full minimal
                )
cpack_add_component(doc
                DISPLAY_NAME "Documentation"
                DESCRIPTION "Class and tool documentation. With tutorials."
                INSTALL_TYPES recommended full
                )
cpack_add_component_group(thirdparty
                     DISPLAY_NAME "Thirdparty binaries"
                     DESCRIPTION "Binaries and files for thirdparty tools and engines."
                     EXPANDED
                     )
foreach(component IN LISTS ${THIRDPARTY_COMPONENT_GROUP})
    cpack_add_component(${component}
                    DISPLAY_NAME ${component}
                    DESCRIPTION "Thirdparty engine ${component}"
                    GROUP thirdparty
                    INSTALL_TYPES recommended full
                    )
endforeach()
