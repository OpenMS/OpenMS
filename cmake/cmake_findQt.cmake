# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Finds the Qt6 modules of the GUI library (WITH_GUI). Included by
# cmake/cmake_findExternalLibs.cmake and, for a build of the GUI against an
# installed OpenMS library, by cmake/find_installed_openms.cmake.
#------------------------------------------------------------------------------
# QT (only needed for GUI)
#------------------------------------------------------------------------------
SET(QT_MIN_VERSION "6.1.0")

if (WITH_GUI)
  find_package(Qt6 ${QT_MIN_VERSION} COMPONENTS Core QUIET)

  IF (Qt6Core_FOUND)
    message(STATUS "Found Qt ${Qt6Core_VERSION}")
  ELSE()
    message(FATAL_ERROR "Qt6Core not found — required when WITH_GUI=ON. Use -DWITH_GUI=OFF to build without GUI.")
  ENDIF()

  # --------------------------------------------------------------------------
  # Find additional Qt libs
  #---------------------------------------------------------------------------
  set (TEMP_OpenMS_GUI_QT_COMPONENTS Gui Widgets Svg OpenGLWidgets)

  # On macOS the platform plugin of QT requires PrintSupport. We link
  # so it's packaged via the bundling/dependency tools/scripts
  if (APPLE)
    set (TEMP_OpenMS_GUI_QT_COMPONENTS ${TEMP_OpenMS_GUI_QT_COMPONENTS} PrintSupport)
  endif()

  set(OpenMS_GUI_QT_COMPONENTS ${TEMP_OpenMS_GUI_QT_COMPONENTS} CACHE INTERNAL "QT components for GUI lib")

  if(NOT NO_WEBENGINE_WIDGETS)
    set(OpenMS_GUI_QT_COMPONENTS_OPT WebEngineWidgets)
  endif()

  find_package(Qt6 REQUIRED COMPONENTS ${OpenMS_GUI_QT_COMPONENTS})

  IF (NOT Qt6Widgets_FOUND OR NOT Qt6Gui_FOUND OR NOT Qt6Svg_FOUND)
    message(STATUS "Qt6Widgets not found!")
    message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -DCMAKE_PREFIX_PATH='<path_to_parent_folder_of_lib_folder_withAllQt6Libs>' <src-dir>")
  ENDIF()
  find_package(Qt6 QUIET COMPONENTS ${OpenMS_GUI_QT_COMPONENTS_OPT})

  # TODO only works if WebEngineWidgets is the only optional component
  set(OpenMS_GUI_QT_FOUND_COMPONENTS_OPT)
  if(Qt6WebEngineWidgets_FOUND)
    list(APPEND OpenMS_GUI_QT_FOUND_COMPONENTS_OPT "WebEngineWidgets")
  else()
    message(WARNING "Qt6WebEngineWidgets not found or disabled, disabling JS Views in TOPPView!")
  endif()

  # The GUI applications derive from TOPPBase and discover tools through ToolHandler,
  # so the tool framework is part of the GUI library's public link interface.
  set(OpenMS_GUI_DEP_LIBRARIES "OpenMS" "OpenMS_CLI")

  foreach(COMP IN LISTS OpenMS_GUI_QT_COMPONENTS)
    list(APPEND OpenMS_GUI_DEP_LIBRARIES "Qt6::${COMP}")
  endforeach()

  foreach(COMP IN LISTS OpenMS_GUI_QT_FOUND_COMPONENTS_OPT)
    list(APPEND OpenMS_GUI_DEP_LIBRARIES "Qt6::${COMP}")
  endforeach()

endif()
