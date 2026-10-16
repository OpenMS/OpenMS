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
## Every TOPP tool links the TOPP tool framework, so the binaries do not start
## without the CLI layer; TOPPView, TOPPAS and the other GUI applications need the
## GUI layer on top of it (library_gui DEPENDS library_cli). Selecting the
## binaries without their libraries installs tools that fail at startup with a
## loader error for libOpenMS_CLI.
## The GUI applications have a component of their own, GUIApplications, except in
## the productbuild installer, which keeps them in Applications
## (OPENMS_GUI_APPLICATIONS_COMPONENT in the top-level CMakeLists.txt).
if(WITH_GUI AND OPENMS_GUI_APPLICATIONS_COMPONENT STREQUAL "Applications")
  set(_openms_applications_depends library_gui)
  set(_openms_applications_description "OpenMS binaries including TOPP tools, TOPPView and TOPPAS.")
else()
  set(_openms_applications_depends library_cli)
  set(_openms_applications_description "The TOPP tools.")
endif()
## The macOS pkg hands pkgbuild a component plist for the app bundles, which keeps the
## installer from relocating them (cmake/generate_applications_component_plist.cmake).
set(_openms_applications_plist)
if(APPLICATIONS_COMPONENT_PLIST)
  set(_openms_applications_plist PLIST "${APPLICATIONS_COMPONENT_PLIST}")
endif()
## Capitalized to match the name install_tool() registers (cmake/install_macros.cmake).
## CPack folds the name to upper case for the CPACK_COMPONENT_<NAME>_* metadata below,
## but compares it verbatim when selecting what to install, so the two must agree.
cpack_add_component(Applications
                DISPLAY_NAME "OpenMS binaries"
                DESCRIPTION "${_openms_applications_description}"
                DEPENDS ${_openms_applications_depends}
                INSTALL_TYPES recommended full minimal
                ${_openms_applications_plist}
                )
unset(_openms_applications_depends)
unset(_openms_applications_description)
unset(_openms_applications_plist)
if(WITH_GUI AND NOT OPENMS_GUI_APPLICATIONS_COMPONENT STREQUAL "Applications")
  cpack_add_component(${OPENMS_GUI_APPLICATIONS_COMPONENT}
                  DISPLAY_NAME "OpenMS GUI applications"
                  DESCRIPTION "TOPPView, TOPPAS and the other GUI applications, and the TOPP tools that need the GUI library."
                  DEPENDS library_gui
                  INSTALL_TYPES recommended full minimal
                  )
endif()
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
