cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

### Installed-consumer test of an OpenMS installation whose prefix was renamed
### after `cmake --install`: what makes OpenMSConfig.cmake *relocatable*
### (OpenMS/OpenMS#4603). Everything the package reports has to be derived from
### the location of the config file itself, because that is what unpacking a
### binary package into an arbitrary directory produces - the prefix a user
### installs to is never the prefix the release was built with.
###
### Driven by the TestExternalCodeRelocated_* tests of src/tests/CMakeLists.txt
### (OPENMS_TEST_INSTALLED_CONSUMER). For an example of an external project see
### ../CMakeLists.txt; a plain consumer project needs none of the checks below.
###
### The installation this is pointed at is the core-only one without the share
### component: OPENMS_DATA_DIR relocates through the same PACKAGE_PREFIX_DIR
### arithmetic as the other directory variables checked here, so copying the
### ~200 MB of shared data a third time would cost more than it proves, and
### ../CMakeLists.txt already checks the installed data of a prefix that stayed
### put. Nothing here may therefore touch <prefix>/share/OpenMS.
project("OpenMS_RelocatedConsumerTest" CXX)

include(CTest)

## The prefix the installation was moved away from, and the one it lives at now.
if(NOT OPENMS_RELOCATED_FROM OR NOT OPENMS_RELOCATED_TO)
  message(FATAL_ERROR "OPENMS_RELOCATED_FROM and OPENMS_RELOCATED_TO have to name the prefix the installation was installed to and the one it was renamed to")
endif()
## Without this the test quietly degenerates into "a second copy of the
## installation also works", which a path baked in at install time would pass:
## it would still resolve, to the prefix that is supposed to be gone.
if(EXISTS "${OPENMS_RELOCATED_FROM}")
  message(FATAL_ERROR "the installation was not moved: '${OPENMS_RELOCATED_FROM}' still exists")
endif()

find_package(OpenMS CONFIG REQUIRED)
message(STATUS "Found OpenMS ${OpenMS_VERSION} at ${OpenMS_DIR}")

#------------------------------------------------------------------------------
## What the package reports has to describe the prefix it was found in.

## Every path of the installation the package names lies below the new prefix.
## This is the check that covers the directories the component set above does not
## install (OPENMS_DATA_DIR, OPENMS_BIN_DIR, OPENMS_DOC_DIR): they cannot be
## tested for existence, but they can be tested for pointing at the right prefix.
## Compared as paths rather than as strings, and safe to compare literally
## because both sides originate from the same string: the driving test passes
## OPENMS_RELOCATED_TO and derives OpenMS_DIR, which @PACKAGE_INIT@ turns back
## into PACKAGE_PREFIX_DIR, from it.
function(openms_assert_below_prefix what path)
  if(path STREQUAL "")
    message(FATAL_ERROR "${what} is empty")
  endif()
  cmake_path(IS_PREFIX OPENMS_RELOCATED_TO "${path}" NORMALIZE _below)
  if(NOT _below)
    message(FATAL_ERROR "${what}='${path}' is not below the prefix the installation was moved to ('${OPENMS_RELOCATED_TO}')")
  endif()
endfunction()

foreach(_dir_var IN ITEMS OPENMS_DATA_DIR OPENMS_SHARE_DIR OPENMS_LIB_DIR OPENMS_BIN_DIR OPENMS_DOC_DIR)
  openms_assert_below_prefix(${_dir_var} "${${_dir_var}}")
endforeach()

## The parts that *are* installed exist at the new location. Since the prefix
## they were installed to is gone (asserted above), existence alone says the
## package resolved them relative to itself.
if(NOT IS_DIRECTORY "${OPENMS_LIB_DIR}")
  message(FATAL_ERROR "OPENMS_LIB_DIR='${OPENMS_LIB_DIR}' is not a directory after the move")
endif()
file(REAL_PATH "${OPENMS_LIB_DIR}" _lib_dir)

foreach(_target IN ITEMS OpenMS::OpenMS OpenMS::OpenSwathAlgo)
  if(NOT TARGET ${_target})
    message(FATAL_ERROR "the relocated installation did not provide the target ${_target}")
  endif()
  ## the imported libraries ...
  get_target_property(_configurations ${_target} IMPORTED_CONFIGURATIONS)
  if(NOT _configurations)
    message(FATAL_ERROR "${_target} has no imported configurations")
  endif()
  foreach(_configuration IN LISTS _configurations)
    get_target_property(_location ${_target} IMPORTED_LOCATION_${_configuration})
    if(NOT EXISTS "${_location}")
      message(FATAL_ERROR "the ${_configuration} location of ${_target} does not exist after the move: ${_location}")
    endif()
    get_filename_component(_location_dir "${_location}" DIRECTORY)
    file(REAL_PATH "${_location_dir}" _location_dir)
    if(NOT _location_dir STREQUAL _lib_dir)
      message(FATAL_ERROR "OPENMS_LIB_DIR='${OPENMS_LIB_DIR}' but the ${_configuration} ${_target} is ${_location}")
    endif()
  endforeach()
  ## ... and the headers they are used with
  get_target_property(_includes ${_target} INTERFACE_INCLUDE_DIRECTORIES)
  if(NOT _includes)
    set(_includes)
  endif()
  foreach(_include IN LISTS _includes)
    openms_assert_below_prefix("an include directory of ${_target}" "${_include}")
    if(NOT IS_DIRECTORY "${_include}")
      message(FATAL_ERROR "the include directory ${_include} of ${_target} does not exist after the move")
    endif()
  endforeach()
endforeach()

#------------------------------------------------------------------------------
## And a program of an external project builds against the moved installation
## and loads its libraries at run time.
add_executable(TestExternalCodeRelocated ../TestExternalCodeRelocated.cpp)
## ${CMAKE_DL_LIBS} for the loader query in the program: dl_iterate_phdr() is part
## of libc on current glibc, but was in libdl before 2.34 and the manylinux images
## of the wheel builds sit right at that boundary. Empty on Windows and macOS.
target_link_libraries(TestExternalCodeRelocated OpenMS::OpenMS ${CMAKE_DL_LIBS})
## What the package reports, for the program to check the library it actually
## loaded against: its version, and the directory the package says its libraries
## are in. The latter is the run-time half of the checks above -- they establish
## what the package *says*, this establishes that the loader agrees.
target_compile_definitions(TestExternalCodeRelocated PRIVATE
                           OPENMS_EXPECTED_VERSION="${OpenMS_VERSION}"
                           OPENMS_EXPECTED_LIB_DIR="${OPENMS_LIB_DIR}"
                           ## disable auto-linking of boost (mainly needed for Windows)
                           BOOST_ALL_NO_LIB)
add_test(TestExternalCodeRelocated TestExternalCodeRelocated)

## Windows has no RPATH: the executable finds OpenMS.dll through PATH.
if(WIN32)
  set_property(TEST TestExternalCodeRelocated APPEND PROPERTY
    ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${OPENMS_LIB_DIR}")
endif()
