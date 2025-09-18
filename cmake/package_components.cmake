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
                DESCRIPTION "Libraries"
                INSTALL_TYPES recommended full minimal
                )
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
