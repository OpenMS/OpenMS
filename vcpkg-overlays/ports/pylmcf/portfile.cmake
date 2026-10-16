vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO michalsta/pylmcf
    REF v0.9.8
    SHA512 db9101ff78f87b797ad8ba2f50432178535e3bc0d03b13d2e987d1d92b01655561930b279ed9a65c365059e94b8a7c2ff3b46fe63b38c8223a28e4335529eeeb
    HEAD_REF main
)

file(INSTALL "${SOURCE_PATH}/src/pylmcf/cpp/"
    DESTINATION "${CURRENT_PACKAGES_DIR}/include")

vcpkg_install_copyright(FILE_LIST
    "${SOURCE_PATH}/LICENCE"
    "${SOURCE_PATH}/src/pylmcf/cpp/lemon/LICENCE"
)

set(VCPKG_POLICY_EMPTY_PACKAGE enabled)
