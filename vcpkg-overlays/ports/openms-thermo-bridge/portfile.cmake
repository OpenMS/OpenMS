if(NOT VCPKG_TARGET_IS_WINDOWS)
    vcpkg_check_linkage(ONLY_DYNAMIC_LIBRARY)
endif()

# Pin both native and managed sources to the same reviewed revision, the one the
# 1.0.0 release (v1.0.0 tag pending) points at, and build the managed component locally.
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO OpenMS/openms-thermo-bridge
    REF e6f0f979732a052ad38d378c2f58790e450770cf
    SHA512 5b7e1b209e1a5efeba90d15a2df99be90995a0630d811e7b9cf4a28ee95c1ba08e18df458cea483748acb395bef7a41df53fddcf1438e047ee9dc95fafbd4865
    HEAD_REF main
    PATCHES
        vcpkg-nethost-use.patch
)

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

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
    -DBUILD_TESTING=OFF
    "-DDOTNET_EXECUTABLE=${DOTNET_EXECUTABLE}"
    -DOPENMS_THERMO_BRIDGE_DOWNLOAD_PREBUILT_MANAGED=OFF
    -DOPENMS_THERMO_BRIDGE_ENABLE_VENDOR_DOWNLOAD=ON
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

# The copyright file above concatenates the notice (a Word document) with the
# bridge's LICENSE, so it cannot be shipped as the notice itself. Install the
# document next to the CMake package as well; OpenMS ships it from there, as it
# does for a bridge built from source (cmake/cmake_findExternalLibs.cmake).
file(INSTALL "${THERMO_LICENSE_PATH}" DESTINATION "${CURRENT_PACKAGES_DIR}/share/OpenMSThermoBridge")
