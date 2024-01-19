# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.12)
if (TARGET tdl)
    return ()
endif ()

option (INSTALL_TDL "Enable installation of TDL. (Projects embedding TDL may want to turn this OFF.)" ON)

find_package (yaml-cpp 0.8.0 QUIET)
if (NOT yaml-cpp_FOUND)
    if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
        message (STATUS "Fetching yaml-cpp")
    endif ()

    include (FetchContent)
    FetchContent_Declare (
        yaml-cpp_fetch_content
        GIT_REPOSITORY "https://github.com/jbeder/yaml-cpp.git"
        # !WORKAROUND Points to first commit after 0.8.0: Fixes CMake deprecation warnings (cmake_minimum_required)
        GIT_TAG "c2680200486572baf8221ba052ef50b58ecd816e")
    option (YAML_CPP_BUILD_CONTRIB "" OFF)
    option (YAML_CPP_BUILD_TOOLS "" OFF)
    option (YAML_BUILD_SHARED_LIBS "" OFF)
    option (YAML_CPP_INSTALL "" ${INSTALL_TDL})
    option (YAML_CPP_BUILD_TESTS "" OFF)
    set_property (GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)
    FetchContent_MakeAvailable (yaml-cpp_fetch_content)
else ()
    if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
        message (STATUS "Found yaml-cpp")
    endif ()
endif ()

add_library (tdl INTERFACE)
target_include_directories (tdl INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src>"
                                          "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries (tdl INTERFACE yaml-cpp::yaml-cpp)
add_library (tdl::tdl ALIAS tdl)

if (INSTALL_TDL)
    include (${CMAKE_CURRENT_LIST_DIR}/cmake/install.cmake)
endif ()
