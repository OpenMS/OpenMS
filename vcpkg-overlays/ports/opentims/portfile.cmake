vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO michalsta/opentims
    REF "v${VERSION}"
    SHA512 43f1b8345a41ac99e4b1ee4c09308f1d9f806bfa91f36c575d52cf7c2757c64877188d65c555c71c357da8dc4efe4b173a81de5b5b0bb0b0ee784b7dd7620a29
    HEAD_REF master
    PATCHES
        fix-static-sqlite-include.patch
)

# Map vcpkg linkage onto upstream's BUILD_SHARED_LIBS switch. Static opentims
# matches OpenMS's FetchContent behavior; shared opentims absorbs sqlite/zstd.
if(VCPKG_LIBRARY_LINKAGE STREQUAL "dynamic")
    set(OPENTIMS_BUILD_SHARED ON)
    set(OPENTIMS_LINK_SQLITE_STATICALLY OFF)
else()
    set(OPENTIMS_BUILD_SHARED OFF)
    set(OPENTIMS_LINK_SQLITE_STATICALLY ON)
endif()

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
    -DBUILD_SHARED_LIBS=${OPENTIMS_BUILD_SHARED}
    -DOPENTIMS_BUILD_LIB=ON
    -DOPENTIMS_BUILD_PYTHON=OFF
    -DOPENTIMS_LINK_SQLITE_STATICALLY=${OPENTIMS_LINK_SQLITE_STATICALLY}
    -DCMAKE_INSTALL_LIBDIR=lib
    -DCMAKE_INSTALL_BINDIR=bin
    -DCMAKE_INSTALL_INCLUDEDIR=include
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME opentims CONFIG_PATH lib/cmake/opentims)
vcpkg_fixup_pkgconfig()

file(REMOVE_RECURSE
    "${CURRENT_PACKAGES_DIR}/debug/include"
    "${CURRENT_PACKAGES_DIR}/share/doc"
    "${CURRENT_PACKAGES_DIR}/debug/share"
)

file(INSTALL "${CMAKE_CURRENT_LIST_DIR}/usage" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENCE")