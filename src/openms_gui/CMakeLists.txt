# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS_GUI")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(CMAKE_AUTOUIC ON) # convert .ui files to ui_....h files using CMake + Qt's UIC
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# --------------------------------------------------------------------------
# include custom macros for openms_gui lib
#include(qt_wrap_ui.cmake)
include(add_mac_bundle.cmake)

# --------------------------------------------------------------------------
# big include file for headers and cpp files
include (${PROJECT_SOURCE_DIR}/includes.cmake)

#------------------------------------------------------------------------------
# add the library
openms_add_library(TARGET_NAME OpenMS_GUI
                   SOURCE_FILES ${OpenMSVisual_sources}
                   HEADER_FILES ${OpenMSVisual_sources_h}
                   INTERNAL_INCLUDES ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include
                   EXTERNAL_INCLUDES
                   LINK_LIBRARIES ${OpenMS_GUI_DEP_LIBRARIES}
                   DLL_EXPORT_PATH "OpenMS/VISUAL/")

# --------------------------------------------------------------------------
# additional linker flags required by openms_gui
if (MSVC)
  ## treat warning of unused function parameter as error, similar to -Werror=unused-variable on GCC
  target_compile_options(OpenMS_GUI PRIVATE "/we4100")
  ## treat warning of unused local variable as error, similar to -Werror=unused-variable on GCC
  target_compile_options(OpenMS_GUI PRIVATE "/we4189")

  set (GUI_lnk_flags "/FORCE:MULTIPLE /INCREMENTAL:NO /ignore:4006 /ignore:4088")
  set_property(TARGET OpenMS_GUI APPEND PROPERTY LINK_FLAGS ${GUI_lnk_flags}) ## allow multiple definitions of symbols (e.g. from template instanciations or STL-derived classes)
endif()


# --------------------------------------------------------------------------
# some targets using OpenMS_GUI (e.g. in the GUI tests) need the ui_..h files, so add them to OpenMS_GUI's interface.
get_property(_isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(_isMultiConfig)
  set(AUTOGEN_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/OpenMS_GUI_autogen/include_$<CONFIG>)
else()
  set(AUTOGEN_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/OpenMS_GUI_autogen/include)
endif()
# add include-directories for UIC'd files (CMake does not do that automatically)
target_include_directories(OpenMS_GUI INTERFACE $<BUILD_INTERFACE:${AUTOGEN_INCLUDE_DIR}>)


#------------------------------------------------------------------------------
# register relevant paths for the doxygen doc generation
openms_doc_path("${PROJECT_SOURCE_DIR}/include")
openms_doc_path("${PROJECT_SOURCE_DIR}/source/VISUAL/APPLICATIONS")

# --------------------------------------------------------------------------
# GUI tools
set(GUI_DIR source/VISUAL/APPLICATIONS/GUITOOLS)
include(${GUI_DIR}/executables.cmake)
foreach(i ${GUI_executables})
  set(resource_file ${PROJECT_SOURCE_DIR}/${GUI_DIR}/${i}.rc)
  set(resource_dir ${PROJECT_SOURCE_DIR}/${GUI_DIR}/${i}-resources/)
  ## add icons to TOPPView, INIFileEditor, etc.
  if (MSVC AND EXISTS ${resource_file})
    message(STATUS "Setting resource file ${resource_file} for ${i}")
    add_executable(${i} ${GUI_DIR}/${i}.cpp ${resource_file})
    install_tool(${i})
  elseif (APPLE AND EXISTS ${resource_dir})
    add_mac_app_bundle(${i})
  else()
    if (MSVC)
      message(STATUS "No resource file (${resource_file}) found for ${i}. No icon will be embedded.")
    elseif (APPLE)
      message(STATUS "No icon file (${i}.icns) found in ${resource_dir}. Will not build ${i} as app bundle.")
    endif()
    add_executable(${i} ${GUI_DIR}/${i}.cpp)
    install_tool(${i})
  endif()

  ## append visual lib as dependency for GUI tools
  target_link_libraries(${i} OpenMS_GUI)
endforeach(i)

  ## install platform plugins. On Mac this is handled for every GUI app inside its app bundle.
  if(WIN32)
      install(IMPORTED_RUNTIME_ARTIFACTS "Qt6::QWindowsIntegrationPlugin"
              DESTINATION "${INSTALL_PLUGIN_DIR}/platforms"
              RUNTIME_DEPENDENCY_SET OPENMS_DEPS
              COMPONENT Dependencies)
      file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
        "[Paths]\nPlugins = ${QT_PLUGIN_PATH_REL_TO_BIN}\n")
      install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
        DESTINATION ${INSTALL_BIN_DIR}
        COMPONENT QtPlatformPlugin)
  endif()

add_custom_target(GUI)
add_dependencies(GUI ${GUI_executables})
