vcpkg_check_linkage(ONLY_STATIC_LIBRARY)

set(EOL_BSPLINE_EXTERN_DIR "${CURRENT_PORT_DIR}/../../../src/openms/extern/eol-bspline")
get_filename_component(EOL_BSPLINE_EXTERN_DIR "${EOL_BSPLINE_EXTERN_DIR}" ABSOLUTE)

if(NOT EXISTS "${EOL_BSPLINE_EXTERN_DIR}/BSpline/BSpline.h")
    message(FATAL_ERROR "eol-bspline sources not found at ${EOL_BSPLINE_EXTERN_DIR}")
endif()

set(SOURCE_PATH "${CURRENT_BUILDTREES_DIR}/src/eol-bspline-staged")
file(REMOVE_RECURSE "${SOURCE_PATH}")
file(MAKE_DIRECTORY "${SOURCE_PATH}")
file(COPY "${EOL_BSPLINE_EXTERN_DIR}/" DESTINATION "${SOURCE_PATH}")

vcpkg_apply_patches(
    SOURCE_PATH "${SOURCE_PATH}"
    PATCHES
        eol-bspline-standalone-cmake.patch
        fix-missing-includes.patch
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
)

vcpkg_cmake_install()   
vcpkg_cmake_config_fixup(PACKAGE_NAME eol-bspline CONFIG_PATH share/eol-bspline)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE" "${SOURCE_PATH}/COPYRIGHT")