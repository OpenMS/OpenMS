if(NOT VCPKG_TARGET_IS_WINDOWS)
    vcpkg_check_linkage(ONLY_DYNAMIC_LIBRARY)
endif()

# Pin the native sources to the v1.1.0 release commit.
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO OpenMS/openms-thermo-bridge
    REF 1d12ef8ab841715f2a223f874d3e9445854ee284
    SHA512 41a2c4ce1977e332d93327f225a443d870012ba0b103d0cd11912f23d9b1db9819da6c292c2c13fd9062b8fec0903b5611efb3d382c8237ec02edd07ebb5fe53
    HEAD_REF main
    PATCHES
        vcpkg-nethost-use.patch
)

# The managed half of the bridge (ThermoWrapperManaged.dll, its runtimeconfig.json and
# the Thermo CommonCore assemblies) comes from the zip that the bridge publishes for each
# platform with every release, verified against the SHA-512 pinned here. Building it in
# this port instead means 'dotnet publish', which downloads five Thermo .nupkg files from
# raw.githubusercontent.com in a single attempt and restores their dependencies from
# nuget.org; that download failed builds at random (OpenMS issue #9717).
#
# OPENMS_THERMO_BRIDGE_RELEASE_TAG is the release tag that points at REF above. When REF
# moves to a commit that has no release, set the tag to "" (and leave the SHA-512 of a
# platform without a published zip empty): the port then publishes the managed half
# itself, which needs a .NET 8 SDK on the build host.
#
# tools/ci/fetch_thermo_assets.sh pins the same zips (as SHA-256) for the builds that do
# not take the bridge from vcpkg; update both together.
set(OPENMS_THERMO_BRIDGE_RELEASE_TAG "v1.1.0")
set(_openms_thermo_bridge_managed_platform "")
set(_openms_thermo_bridge_managed_sha512 "")
if(VCPKG_TARGET_IS_LINUX AND VCPKG_TARGET_ARCHITECTURE STREQUAL "x64")
    set(_openms_thermo_bridge_managed_platform "linux-x64")
    set(_openms_thermo_bridge_managed_sha512 "e1baba71106a988f5f5d2b9a21df4ab9910292331f57579562a742d84dfd384ac48bc0c8b7a7f191b123156b27c1796a141e784c2ebd85b3a811982f6ceee41b")
elseif(VCPKG_TARGET_IS_LINUX AND VCPKG_TARGET_ARCHITECTURE STREQUAL "arm64")
    set(_openms_thermo_bridge_managed_platform "linux-arm64")
    set(_openms_thermo_bridge_managed_sha512 "85209470986c46cfaf43aa11d085dd3b352fb7946ace7ad2b5b0d87363ce1970a03e4da65eb8839e9a5bc196fd5c2e9df5c5b1e5296f4a0e3f8b420e52ae4d7b")
elseif(VCPKG_TARGET_IS_OSX AND VCPKG_TARGET_ARCHITECTURE STREQUAL "arm64")
    set(_openms_thermo_bridge_managed_platform "osx-arm64")
    set(_openms_thermo_bridge_managed_sha512 "d6dbd695c576e0daf8552dea7f5bf8f92e6ca9d7d99c8351ce636ca2852149719b48d3bd30e6c730de195dc05e26250d59fefb5cb578fb89f286902cf0f46eb7")
elseif(VCPKG_TARGET_IS_WINDOWS AND VCPKG_TARGET_ARCHITECTURE STREQUAL "x64")
    set(_openms_thermo_bridge_managed_platform "win-x64")
    set(_openms_thermo_bridge_managed_sha512 "71e26c309138c0e2d3afd754dc0b65a82402d9f48d0b7aed168d61873fabb053d80841bef8e9f4f8a296aa3ad520e92e068755d9c967399197dde59cc3393b04")
endif()

if(OPENMS_THERMO_BRIDGE_RELEASE_TAG AND _openms_thermo_bridge_managed_sha512)
    set(_openms_thermo_bridge_managed_zip
        "openms-thermo-bridge-managed-${_openms_thermo_bridge_managed_platform}-${OPENMS_THERMO_BRIDGE_RELEASE_TAG}.zip")
    vcpkg_download_distfile(_openms_thermo_bridge_managed_archive
        URLS "https://github.com/OpenMS/openms-thermo-bridge/releases/download/${OPENMS_THERMO_BRIDGE_RELEASE_TAG}/${_openms_thermo_bridge_managed_zip}"
        FILENAME "${_openms_thermo_bridge_managed_zip}"
        SHA512 "${_openms_thermo_bridge_managed_sha512}"
    )
    # The zip holds the assemblies at its top level, so there is no directory to strip.
    vcpkg_extract_source_archive(_openms_thermo_bridge_managed_dir
        ARCHIVE "${_openms_thermo_bridge_managed_archive}"
        SOURCE_BASE managed
        NO_REMOVE_ONE_LEVEL
    )
    set(_openms_thermo_bridge_managed_options
        "-DOPENMS_THERMO_BRIDGE_PREBUILT_MANAGED_DIR=${_openms_thermo_bridge_managed_dir}")
else()
    message(STATUS "openms-thermo-bridge: no pre-built managed bridge is pinned for "
                   "${TARGET_TRIPLET}; publishing it with the .NET SDK")

    find_program(DOTNET_EXECUTABLE NAMES dotnet
        HINTS "$ENV{DOTNET_ROOT}" "$ENV{ProgramFiles}/dotnet"
              "/opt/homebrew/opt/dotnet/bin" REQUIRED)

    # ThermoWrapperManaged.csproj targets net8.0, so 'dotnet publish' needs an SDK of
    # at least that major version; the bridge itself only checks that some SDK exists.
    execute_process(
        COMMAND "${DOTNET_EXECUTABLE}" --list-sdks
        OUTPUT_VARIABLE _openms_thermo_bridge_sdks
        ERROR_VARIABLE _openms_thermo_bridge_sdks_error
        RESULT_VARIABLE _openms_thermo_bridge_sdks_result
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(_openms_thermo_bridge_has_net8 FALSE)
    if(_openms_thermo_bridge_sdks_result EQUAL 0)
        string(REPLACE "\n" ";" _openms_thermo_bridge_sdk_lines "${_openms_thermo_bridge_sdks}")
        foreach(_line IN LISTS _openms_thermo_bridge_sdk_lines)
            # each line reads "<major>.<minor>.<patch> [<path>]"
            if(_line MATCHES "^([0-9]+)\\." AND CMAKE_MATCH_1 GREATER_EQUAL 8)
                set(_openms_thermo_bridge_has_net8 TRUE)
            endif()
        endforeach()
    endif()
    if(NOT _openms_thermo_bridge_has_net8)
        message(FATAL_ERROR
            "openms-thermo-bridge: building the managed bridge requires a .NET SDK 8.0 or newer "
            "('${DOTNET_EXECUTABLE} --list-sdks' reported: '${_openms_thermo_bridge_sdks}' "
            "${_openms_thermo_bridge_sdks_error}). Install the .NET 8 SDK or set DOTNET_ROOT.")
    endif()

    set(_openms_thermo_bridge_managed_options
        "-DDOTNET_EXECUTABLE=${DOTNET_EXECUTABLE}"
        -DOPENMS_THERMO_BRIDGE_ENABLE_VENDOR_DOWNLOAD=ON)
endif()

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
    -DBUILD_TESTING=OFF
    -DOPENMS_THERMO_BRIDGE_DOWNLOAD_PREBUILT_MANAGED=OFF
    -DOPENMS_THERMO_BRIDGE_BUILD_CLI=OFF
    ${_openms_thermo_bridge_managed_options}
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

# The copyright file above concatenates the notice (a Word document) with the
# bridge's LICENSE, so it cannot be shipped as the notice itself. Install the
# document next to the CMake package as well; OpenMS ships it from there, as it
# does for a bridge built from source (cmake/cmake_findExternalLibs.cmake).
file(INSTALL "${THERMO_LICENSE_PATH}" DESTINATION "${CURRENT_PACKAGES_DIR}/share/OpenMSThermoBridge")
