if(NOT VCPKG_TARGET_IS_WINDOWS)
    vcpkg_check_linkage(ONLY_DYNAMIC_LIBRARY)
endif()

vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO OpenMS/openms-thermo-bridge
    REF "v${VERSION}"
    SHA512 436daff0a332a48da3945eb85e03a49612a0b2489dbbf64692c5ff1b8ac6df00c1cad2fe9fd987166662ab675e8a234c1a026911c47abac9ca6b8f9cf7cb7151
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

file(DOWNLOAD
    "https://raw.githubusercontent.com/thermofisherlsms/RawFileReader/80963674b5c10e58236da63023ad6fa0264bbb00/License.doc"
    "${CURRENT_BUILDTREES_DIR}/share-tmp/ThermoRawFileReader-License.doc"
)
vcpkg_install_copyright(FILE_LIST "${CURRENT_BUILDTREES_DIR}/share-tmp/ThermoRawFileReader-License.doc")