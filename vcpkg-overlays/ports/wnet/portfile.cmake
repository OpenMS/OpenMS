vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO michalsta/wnet
    REF v0.9.11
    SHA512 6ea620ddde04402df1112396ea217576ba01858d693ace02980a2555b9748f1d92d4cc9a75fea6642e55390af004f57c1f2ced8bb1f97f8b3a3d2957d7c836ba
    HEAD_REF main
)

file(INSTALL "${SOURCE_PATH}/src/wnet/cpp/"
    DESTINATION "${CURRENT_PACKAGES_DIR}/include")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENCE")

set(VCPKG_POLICY_EMPTY_PACKAGE enabled)
