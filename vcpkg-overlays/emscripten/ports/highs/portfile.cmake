vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO ERGO-Code/HiGHS
    REF "v1.14.0"
    SHA512 d4c2d09bf41384a58a04528ffaac223eee40dea4a1cadaea369791f6e8d1aa61ff3fce57e7c97cb5f4bef0f1305681066e8ef7483f4471770295f95b9ec337f7
    HEAD_REF master
)

# OpenMS: for WebAssembly, only the library. The command-line solver would be built as
# highs.js, which vcpkg_copy_tools() does not find; and without threads, HiGHS has to run
# its task scheduler with one worker (which it would also choose for the single core
# std::thread::hardware_concurrency() reports there).
set(OPENMS_HIGHS_OPTIONS "")
if(VCPKG_TARGET_IS_EMSCRIPTEN)
    list(APPEND OPENMS_HIGHS_OPTIONS -DBUILD_CXX_EXE=OFF -DHIGHS_NO_DEFAULT_THREADS=ON)
endif()

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DFAST_BUILD=ON
        -DBUILD_TESTING=OFF
        -DBUILD_EXAMPLES=OFF
        -DCMAKE_REQUIRE_FIND_PACKAGE_ZLIB=ON
        ${OPENMS_HIGHS_OPTIONS}
)

vcpkg_cmake_install()
vcpkg_fixup_pkgconfig()
if(NOT VCPKG_TARGET_IS_EMSCRIPTEN)
    vcpkg_copy_tools(TOOL_NAMES highs AUTO_CLEAN)
endif()

vcpkg_cmake_config_fixup(CONFIG_PATH "lib/cmake/highs")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE.txt")

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share")
