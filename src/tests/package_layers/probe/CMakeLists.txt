cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

### Probe of the package-layer tests: configured against the build tree of the
### stub project in .. and checks what its OpenMSConfig.cmake reports.
project("OpenMS_PackageLayersProbe" CXX)

option(EXPECT_GUI "The build tree was configured with WITH_GUI=ON, so the GUI layer has to be there" OFF)
option(REQUIRE_GUI "Require the GUI component, which a WITH_GUI=OFF build tree has to refuse" OFF)

if(REQUIRE_GUI)
  ## has to fail: the driving test passes on the reason OpenMSConfig.cmake gives
  find_package(OpenMS CONFIG REQUIRED COMPONENTS GUI)
  message(FATAL_ERROR "find_package(OpenMS COMPONENTS GUI) succeeded on a build tree configured with WITH_GUI=OFF")
endif()

## GUI as an optional component: its Qt6 lookup must not decide these checks
find_package(OpenMS CONFIG REQUIRED OPTIONAL_COMPONENTS GUI)

foreach(_target IN ITEMS OpenMS::OpenMS OpenMS::OpenSwathAlgo OpenMS::OpenMS_CLI)
  if(NOT TARGET ${_target})
    message(FATAL_ERROR "the build tree's package did not provide the target ${_target}")
  endif()
endforeach()
if(NOT OpenMS_CLI_FOUND)
  message(FATAL_ERROR "the build tree's package does not report its CLI layer")
endif()

if(EXPECT_GUI)
  if(NOT OpenMS_WITH_GUI OR NOT TARGET OpenMS::OpenMS_GUI)
    message(FATAL_ERROR "the build tree was configured with WITH_GUI=ON but its package does not provide the GUI layer")
  endif()
else()
  foreach(_var IN ITEMS OpenMS_WITH_GUI OpenMS_GUI_FOUND)
    if(${_var})
      message(FATAL_ERROR "the build tree was configured with WITH_GUI=OFF but its package reports ${_var}=${${_var}}")
    endif()
  endforeach()
  foreach(_target IN ITEMS OpenMS::OpenMS_GUI OpenMS_GUI)
    if(TARGET ${_target})
      message(FATAL_ERROR "the build tree was configured with WITH_GUI=OFF but its package provides the target ${_target} (stale GUI target file)")
    endif()
  endforeach()
endif()
message(STATUS "package layers as expected (GUI: ${OpenMS_WITH_GUI})")
