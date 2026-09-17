cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

### Stand-in for the OpenMS build in the package-layer tests of
### src/tests/CMakeLists.txt (OPENMS_TEST_INSTALLED_CONSUMER). It uses the real
### export and install macros and the real package template with stub libraries,
### so the CMake package of its build tree can be checked in seconds, including
### after reconfiguring the same build directory with another WITH_GUI setting.
### Nothing is compiled: the build-tree target files are written at generation.
project("OpenMS_PackageLayersStub" CXX)

option(WITH_GUI "Build the stub GUI library (the GUI layer of the package)" OFF)

# what the macros and cmake/OpenMSConfig.cmake.in expect from the OpenMS build
get_filename_component(OPENMS_HOST_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(OPENMS_HOST_BINARY_DIRECTORY "${PROJECT_BINARY_DIR}")
set(OPENMS_PACKAGE_VERSION 3.6.0)
set(INSTALL_LIB_DIR lib)
set(INSTALL_BIN_DIR bin)
set(INSTALL_INCLUDE_DIR include)
set(INSTALL_SHARE_DIR share/OpenMS)
set(INSTALL_DOC_DIR share/doc)
set(INSTALL_CMAKE_DIR lib/cmake/OpenMS)

include(${OPENMS_HOST_DIRECTORY}/cmake/install_macros.cmake)
include(${OPENMS_HOST_DIRECTORY}/cmake/export_macros.cmake)

# the libraries of the three layers, registered the way openms_add_library() does it
add_library(OpenSwathAlgo SHARED stub.cpp)
add_library(OpenMS SHARED stub.cpp)
target_link_libraries(OpenMS PUBLIC OpenSwathAlgo)
add_library(OpenMS_CLI SHARED stub.cpp)
target_link_libraries(OpenMS_CLI PUBLIC OpenMS)
foreach(_core IN ITEMS OpenSwathAlgo OpenMS)
  install_library(${_core})
  openms_register_export_target(${_core})
endforeach()
install_library(OpenMS_CLI EXPORT_SET ${OPENMS_CLI_EXPORT_SET})
openms_register_export_target(OpenMS_CLI ${OPENMS_CLI_EXPORT_SET})
if(WITH_GUI)
  add_library(OpenMS_GUI SHARED stub.cpp)
  target_link_libraries(OpenMS_GUI PUBLIC OpenMS_CLI)
  install_library(OpenMS_GUI EXPORT_SET ${OPENMS_GUI_EXPORT_SET})
  openms_register_export_target(OpenMS_GUI ${OPENMS_GUI_EXPORT_SET})
endif()
# install_library() names this set; it is never installed here
install(RUNTIME_DEPENDENCY_SET OPENMS_DEPS DESTINATION ${INSTALL_LIB_DIR})

openms_export_targets()
install_export_targets()
