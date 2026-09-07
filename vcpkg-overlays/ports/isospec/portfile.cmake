vcpkg_check_linkage(ONLY_STATIC_LIBRARY)

set(ISOSPEC_EXTERN_DIR "${CURRENT_PORT_DIR}/../../../src/openms/extern/IsoSpec")
get_filename_component(ISOSPEC_EXTERN_DIR "${ISOSPEC_EXTERN_DIR}" ABSOLUTE)

if(NOT EXISTS "${ISOSPEC_EXTERN_DIR}/IsoSpec/isoSpec++.h")
    message(FATAL_ERROR "IsoSpec sources not found at ${ISOSPEC_EXTERN_DIR}")
endif()

set(SOURCE_PATH "${CURRENT_BUILDTREES_DIR}/src/isospec-staged")
file(REMOVE_RECURSE "${SOURCE_PATH}")
file(MAKE_DIRECTORY "${SOURCE_PATH}")
file(COPY "${ISOSPEC_EXTERN_DIR}/" DESTINATION "${SOURCE_PATH}")

vcpkg_apply_patches(
    SOURCE_PATH "${SOURCE_PATH}"
    PATCHES
        isospec-standalone-cmake.patch
        isospec-building-openms-value.patch
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME isospec CONFIG_PATH share/isospec)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

vcpkg_install_copyright(
    FILE_LIST "${SOURCE_PATH}/LICENSE"
    COMMENT
        "IsoSpec itself is licenses under BSD-2-Clause (see LICENSE)
        Additionally, this library bundles some files with their own licenses:
        - mman.h: MIT License. Comes from the repository (https://github.com/witwall/mman-win32). No explicit license file is included, but stated in file header's comment.
        - btrd.h: Boost Software License. Adapted from boost random/binomial_distribution.hpp header file, at version 1.71. Copyright Steven Watanabe 2010. Full license text: http://www.boost.org/LICENSE_1_0.txt"
)