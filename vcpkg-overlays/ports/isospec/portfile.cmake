set(SOURCE_PATH "${CURRENT_PORT_DIR}/src")

vcpkg_check_linkage(ONLY_STATIC_LIBRARY)
vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME isospec CONFIG_PATH share/isospec)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

vcpkg_install_copyright(
    FILE_LIST "${SOURCE_PATH}/IsoSpec/LICENSE"
    COMMENT
        "IsoSpec itself is licenses under BSD-2-Clause (see LICENSE)

        Additionally, this library bundles some files with their own licenses:

        - mman.h: MIT License. Comes from the repository (https://github.com/witwall/mman-win32). No explicit license file is included, but stated in file header's comment.

        - btrd.h: Boost Software License. Adapted from boost random/binomial_distribution.hpp header file, at version 1.71. Copyright Steven Watanabe 2010. Full license text: http://www.boost.org/LICENSE_1_0.txt"
)