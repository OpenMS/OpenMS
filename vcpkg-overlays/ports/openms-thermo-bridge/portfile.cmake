if(NOT VCPKG_TARGET_IS_WINDOWS)
    vcpkg_check_linkage(ONLY_DYNAMIC_LIBRARY)
endif()

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO OpenMS/openms-thermo-bridge
    REF e05aa5853474943a43819151e403325860f00168
    SHA512 b0e9f78dc90d748edba51145b2b1a07d3adf9cb4f3f4274388b008f72b8b927a36ab6f0e10f80155ba58dfdcc3bb405b5f37643af1c92a73dbe58a2db4441fff
    HEAD_REF main
    PATCHES
        vcpkg-nethost-use.patch
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
    -DBUILD_TESTING=OFF
    -DOPENMS_THERMO_BRIDGE_DOWNLOAD_PREBUILT_MANAGED=ON
    -DOPENMS_THERMO_BRIDGE_BUILD_CLI=OFF
)

vcpkg_cmake_install()

vcpkg_cmake_config_fixup(
    PACKAGE_NAME OpenMSThermoBridge
    CONFIG_PATH lib/cmake/OpenMSThermoBridge
)

vcpkg_copy_pdbs()

file(REMOVE_RECURSE
    "${CURRENT_PACKAGES_DIR}/debug/include"
    "${CURRENT_PACKAGES_DIR}/debug/lib/openms_thermo_bridge/managed"
)

file(MAKE_DIRECTORY
    "${CURRENT_BUILDTREES_DIR}/share-tmp"
)

vcpkg_download_distfile(THERMO_LICENSE_PATH
    URLS "https://raw.githubusercontent.com/thermofisherlsms/RawFileReader/80963674b5c10e58236da63023ad6fa0264bbb00/License.doc"
    FILENAME "ThermoRawFileReader-License.doc"
    SHA512 6ecc1691854ebd16914b2035c20585f8d5afd5f6fcf0a4ef3564ee7d6bfe65c8d332df2e56fa926cf093602a7983ae6a5c300f7ac23e3241d9ea1dc9f3b01b03
)

vcpkg_install_copyright(FILE_LIST "${THERMO_LICENSE_PATH}" "${SOURCE_PATH}/LICENSE")