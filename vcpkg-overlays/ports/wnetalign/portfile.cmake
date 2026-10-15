vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO michalsta/wnetalign
    REF v0.9.8
    SHA512 9d6fb75f720e386ff06c31bea89f8f844034268c5a9a6eb4d1d9deeea46d0710863588f354e5f9af38cb2184e9a3feeec110c24f227c4493ec0046557146b235
    HEAD_REF main
)

file(INSTALL "${SOURCE_PATH}/src/wnetalign/cpp/"
    DESTINATION "${CURRENT_PACKAGES_DIR}/include")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENCE")

set(VCPKG_POLICY_EMPTY_PACKAGE enabled)
