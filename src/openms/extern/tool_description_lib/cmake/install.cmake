# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.12)

include (${CMAKE_CURRENT_LIST_DIR}/version.cmake)
include (GNUInstallDirs)

# Needed such that the generated tdl-config.cmake properly links against the also-to-be-installed yaml-cpp.
if (NOT yaml-cpp_FOUND)
    set (TDL_EXPORT_YAML_CPP "yaml-cpp")
endif ()
install (TARGETS tdl ${TDL_EXPORT_YAML_CPP}
         EXPORT tdl_targets
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         FRAMEWORK DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../src/tdl" TYPE INCLUDE)
install (EXPORT tdl_targets
         NAMESPACE tdl::
         DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tdl
         FILE tdl-config.cmake)

include (CMakePackageConfigHelpers)
set (version_file "${CMAKE_CURRENT_BINARY_DIR}/cmake/tdl-config-version.cmake")
write_basic_package_version_file (
    ${version_file}
    VERSION "${TDL_VERSION}"
    COMPATIBILITY AnyNewerVersion)
install (FILES ${version_file} DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tdl)

install (FILES "${CMAKE_CURRENT_LIST_DIR}/../README.md" TYPE DOC)
install (DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../LICENSES" TYPE DOC)
