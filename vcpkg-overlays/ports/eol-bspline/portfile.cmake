set(SOURCE_PATH "${CURRENT_PORT_DIR}/src")

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
)

vcpkg_cmake_install()   
vcpkg_cmake_config_fixup(PACKAGE_NAME eol-bspline CONFIG_PATH share/eol-bspline)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE" "${SOURCE_PATH}/COPYRIGHT")