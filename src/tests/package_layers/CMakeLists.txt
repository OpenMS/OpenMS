cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

### Stand-in for the OpenMS build in the package-layer tests of
### src/tests/CMakeLists.txt (OPENMS_TEST_INSTALLED_CONSUMER). It uses the real
### export and install macros and the real package template with stub libraries,
### so the CMake package of its build tree can be checked in seconds, including
### after reconfiguring the same build directory with another WITH_CLI setting.
### Nothing is compiled: the build-tree target files are written at generation.
###
### With WITH_GUI_PACKAGE it also writes a stub of the OpenMSGUI package of
### OpenMS-TOPP-GUI next to OpenMSConfig.cmake (where the GUI layer puts it on
### Windows), to which OpenMSConfig.cmake delegates its component GUI.
project("OpenMS_PackageLayersStub" CXX)

option(WITH_CLI "Build the stub CLI library (the CLI layer of the package)" ON)
option(WITH_GUI_PACKAGE "Write a stub OpenMSGUI package next to the package of this build tree" OFF)

# what the macros and cmake/OpenMSConfig.cmake.in expect from the OpenMS build
get_filename_component(OPENMS_HOST_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(OPENMS_HOST_BINARY_DIRECTORY "${PROJECT_BINARY_DIR}")
set(OPENMS_PACKAGE_VERSION 3.6.0)
set(OPENMS_PACKAGE_VERSION_FULLSTRING ${OPENMS_PACKAGE_VERSION})
set(OPENMS_PACKAGE_VERSION_PRERELEASE_IDENTIFIER "")
set(INSTALL_LIB_DIR lib)
set(INSTALL_BIN_DIR bin)
set(INSTALL_INCLUDE_DIR include)
set(INSTALL_SHARE_DIR share/OpenMS)
set(INSTALL_DOC_DIR share/doc)
set(INSTALL_CMAKE_DIR lib/cmake/OpenMS)

include(${OPENMS_HOST_DIRECTORY}/cmake/install_macros.cmake)
include(${OPENMS_HOST_DIRECTORY}/cmake/export_macros.cmake)

# the libraries of the two layers, registered the way openms_add_library() does it
add_library(OpenSwathAlgo SHARED stub.cpp)
add_library(OpenMS SHARED stub.cpp)
target_link_libraries(OpenMS PUBLIC OpenSwathAlgo)
foreach(_core IN ITEMS OpenSwathAlgo OpenMS)
  install_library(${_core})
  openms_register_export_target(${_core})
endforeach()
if(WITH_CLI)
  add_library(OpenMS_CLI SHARED stub.cpp)
  target_link_libraries(OpenMS_CLI PUBLIC OpenMS)
  install_library(OpenMS_CLI EXPORT_SET ${OPENMS_CLI_EXPORT_SET})
  openms_register_export_target(OpenMS_CLI ${OPENMS_CLI_EXPORT_SET})
endif()
# install_library() names this set; it is never installed here
install(RUNTIME_DEPENDENCY_SET OPENMS_DEPS DESTINATION ${INSTALL_LIB_DIR})

openms_export_targets()
install_export_targets()

# The stub GUI package: an interface target on top of the CLI layer. Unlike the real
# OpenMSGUIConfig.cmake it refuses to search OpenMS itself, so the probe also checks
# that OpenMSConfig.cmake loads the core and CLI targets before it delegates.
set(_gui_config "${PROJECT_BINARY_DIR}/OpenMSGUIConfig.cmake")
set(_gui_version "${PROJECT_BINARY_DIR}/OpenMSGUIConfigVersion.cmake")
if(WITH_GUI_PACKAGE)
  file(WRITE "${_gui_config}" [=[
if(NOT TARGET OpenMS::OpenMS_CLI)
  set(OpenMSGUI_FOUND FALSE)
  set(OpenMSGUI_NOT_FOUND_MESSAGE "the stub GUI package was loaded before the CLI layer")
  return()
endif()
if(NOT TARGET OpenMS::OpenMS_GUI)
  add_library(OpenMS::OpenMS_GUI INTERFACE IMPORTED)
  set_target_properties(OpenMS::OpenMS_GUI PROPERTIES INTERFACE_LINK_LIBRARIES OpenMS::OpenMS_CLI)
endif()
if(NOT TARGET OpenMS_GUI)
  add_library(OpenMS_GUI ALIAS OpenMS::OpenMS_GUI)
endif()
]=])
  include(CMakePackageConfigHelpers)
  write_basic_package_version_file("${_gui_version}"
    VERSION ${OPENMS_PACKAGE_VERSION} COMPATIBILITY SameMinorVersion ARCH_INDEPENDENT)
else()
  file(REMOVE "${_gui_config}" "${_gui_version}")
endif()
