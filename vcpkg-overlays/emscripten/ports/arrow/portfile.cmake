vcpkg_download_distfile(
    ARCHIVE_PATH
    URLS "https://archive.apache.org/dist/arrow/arrow-${VERSION}/apache-arrow-${VERSION}.tar.gz"
    FILENAME apache-arrow-${VERSION}.tar.gz
    SHA512 ffbac0c198c74f26b2ac173aa2b054cfd3902927ca3db36fbbae25353e8696c32b744e805c25801ced22641af71a5a8ea441202d8d3494aa803b797995d9d408
)
vcpkg_extract_source_archive(
    SOURCE_PATH
    ARCHIVE ${ARCHIVE_PATH}
    PATCHES
        0001-msvc-static-name.patch
        0003-android-musl.patch
        0004-android-datetime.patch
        0005-cmake-msvcruntime.patch
        0007-use-vcpkg-mimalloc.patch
)

# Check cpp/cmake_modules/DefineOptions.cmake for option dependencies -
# they must be modeled as feature dependencies in vcpkg.json.
vcpkg_check_features(OUT_FEATURE_OPTIONS FEATURE_OPTIONS
    FEATURES
        acero       ARROW_ACERO
        compute     ARROW_COMPUTE
        csv         ARROW_CSV
        cuda        ARROW_CUDA
        dataset     ARROW_DATASET
        filesystem  ARROW_FILESYSTEM
        flight      ARROW_FLIGHT
        flightsql   ARROW_FLIGHT_SQL
        gcs         ARROW_GCS
        jemalloc    ARROW_JEMALLOC
        json        ARROW_JSON
        mimalloc    ARROW_MIMALLOC
        orc         ARROW_ORC
        parquet     ARROW_PARQUET
        parquet     PARQUET_REQUIRE_ENCRYPTION
        s3          ARROW_S3
)

if(VCPKG_TARGET_IS_WINDOWS AND NOT VCPKG_TARGET_IS_MINGW)
    list(APPEND FEATURE_OPTIONS "-DARROW_USE_NATIVE_INT128=OFF")
endif()

if(VCPKG_TARGET_IS_WINDOWS AND VCPKG_TARGET_ARCHITECTURE STREQUAL "arm64")
    list(APPEND FEATURE_OPTIONS "-DARROW_SIMD_LEVEL=NONE")
endif()

# OpenMS: Arrow's build adds -fexceptions (exceptions emulated in JavaScript) for
# Emscripten. All code linked into one WebAssembly module must use the same exception
# handling, and the wasm32-emscripten triplet compiles everything with native WebAssembly
# exceptions (-fwasm-exceptions): leave the choice to the triplet's flags.
if(VCPKG_TARGET_IS_EMSCRIPTEN)
    set(_arrow_cxx_flags_file "${SOURCE_PATH}/cpp/cmake_modules/SetupCxxFlags.cmake")
    vcpkg_replace_string("${_arrow_cxx_flags_file}" [[set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -fexceptions")]] [[set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")]])
    vcpkg_replace_string("${_arrow_cxx_flags_file}" [["${CMAKE_CXX_FLAGS} -fPIC -fexceptions -Wno-error=deprecated-literal-operator")]] [["${CMAKE_CXX_FLAGS} -fPIC -Wno-error=deprecated-literal-operator")]])
    vcpkg_replace_string("${_arrow_cxx_flags_file}" [["-sWASM_BIGINT=1 -fexceptions -Wno-error=linkflags"]] [["-sWASM_BIGINT=1 -Wno-error=linkflags"]])
endif()

# OpenMS: the settings of Arrow's own Emscripten preset (ninja-release-emscripten in
# cpp/CMakePresets.json, which Pyodide's pyarrow is built from). A WebAssembly module
# without pthreads cannot start threads, so Arrow's thread pools are compiled out.
# Parquet modular encryption needs OpenSSL, which is not built for WebAssembly; the
# options after the feature options override PARQUET_REQUIRE_ENCRYPTION=ON.
if(VCPKG_TARGET_IS_EMSCRIPTEN)
    list(APPEND FEATURE_OPTIONS
        "-DARROW_ENABLE_THREADING=OFF"
        "-DARROW_SIMD_LEVEL=NONE"
        "-DARROW_RUNTIME_SIMD_LEVEL=NONE"
        "-DARROW_JEMALLOC=OFF"
        "-DARROW_MIMALLOC=OFF"
        "-DPARQUET_REQUIRE_ENCRYPTION=OFF"
    )
endif()

string(COMPARE EQUAL ${VCPKG_LIBRARY_LINKAGE} "dynamic" ARROW_BUILD_SHARED)
string(COMPARE EQUAL ${VCPKG_LIBRARY_LINKAGE} "static" ARROW_BUILD_STATIC)
string(COMPARE EQUAL ${VCPKG_LIBRARY_LINKAGE} "dynamic" ARROW_DEPENDENCY_USE_SHARED)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}/cpp"
    OPTIONS
        ${FEATURE_OPTIONS}
        -DARROW_BUILD_SHARED=${ARROW_BUILD_SHARED}
        -DARROW_BUILD_STATIC=${ARROW_BUILD_STATIC}
        -DARROW_BUILD_TESTS=OFF
        -DARROW_DEPENDENCY_SOURCE=SYSTEM
        -DARROW_DEPENDENCY_USE_SHARED=${ARROW_DEPENDENCY_USE_SHARED}
        -DARROW_PACKAGE_KIND=vcpkg
        -DARROW_WITH_BROTLI=ON
        -DARROW_WITH_BZ2=ON
        -DARROW_WITH_LZ4=ON
        -DARROW_WITH_SNAPPY=ON
        -DARROW_WITH_UTF8PROC=ON
        -DARROW_WITH_ZLIB=ON
        -DARROW_WITH_ZSTD=ON
        -DBUILD_WARNING_LEVEL=PRODUCTION
        -DZSTD_MSVC_LIB_PREFIX=
    MAYBE_UNUSED_VARIABLES
        ZSTD_MSVC_LIB_PREFIX
)

vcpkg_cmake_install()
vcpkg_copy_pdbs()

vcpkg_fixup_pkgconfig()

if(EXISTS "${CURRENT_PACKAGES_DIR}/lib/arrow_static.lib")
    message(FATAL_ERROR "Installed lib file should be named 'arrow.lib' via patching the upstream build.")
endif()

if("dataset" IN_LIST FEATURES)
    vcpkg_cmake_config_fixup(
        PACKAGE_NAME arrowdataset
        CONFIG_PATH lib/cmake/ArrowDataset
        DO_NOT_DELETE_PARENT_CONFIG_PATH
    )
endif()

if("acero" IN_LIST FEATURES)
    vcpkg_cmake_config_fixup(
        PACKAGE_NAME arrowacero
        CONFIG_PATH lib/cmake/ArrowAcero
        DO_NOT_DELETE_PARENT_CONFIG_PATH
    )
endif()

if("compute" IN_LIST FEATURES)
    vcpkg_cmake_config_fixup(
        PACKAGE_NAME arrowcompute
        CONFIG_PATH lib/cmake/ArrowCompute
        DO_NOT_DELETE_PARENT_CONFIG_PATH
    )
endif()

if("flight" IN_LIST FEATURES)
    vcpkg_cmake_config_fixup(
        PACKAGE_NAME arrowflight
        CONFIG_PATH lib/cmake/ArrowFlight
        DO_NOT_DELETE_PARENT_CONFIG_PATH
    )
endif()

if("flightsql" IN_LIST FEATURES)
    vcpkg_cmake_config_fixup(
        PACKAGE_NAME arrowflightsql
        CONFIG_PATH lib/cmake/ArrowFlightSql
        DO_NOT_DELETE_PARENT_CONFIG_PATH
    )
endif()

if("parquet" IN_LIST FEATURES)
    vcpkg_cmake_config_fixup(
        PACKAGE_NAME parquet
        CONFIG_PATH lib/cmake/Parquet
        DO_NOT_DELETE_PARENT_CONFIG_PATH
    )
endif()

file(GLOB main_configs "${CURRENT_PACKAGES_DIR}/lib/cmake/Arrow/*onfig.cmake")
file(GLOB extra_configs "${CURRENT_PACKAGES_DIR}/lib/cmake/*/*onfig.cmake")
list(REMOVE_ITEM extra_configs ${main_configs})
if(NOT "${extra_configs}" STREQUAL "")
    message("${Z_VCPKG_BACKCOMPAT_MESSAGE_LEVEL}"
        "Unhandled CMake config: ${extra_configs}\n"
        "This might be caused by insufficient feature dependencies in ports/arrow/vcpkg.json."
    )
endif()
vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/Arrow)

file(INSTALL "${CMAKE_CURRENT_LIST_DIR}/usage" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}")
foreach(feature IN ITEMS parquet dataset acero compute flight flightsql)
    if(feature IN_LIST FEATURES)
        file(READ "${CMAKE_CURRENT_LIST_DIR}/usage-${feature}" feature_usage)
        file(APPEND "${CURRENT_PACKAGES_DIR}/share/${PORT}/usage" "${feature_usage}")
    endif()
endforeach()

if("example" IN_LIST FEATURES)
    file(INSTALL "${SOURCE_PATH}/cpp/examples/minimal_build/" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}/example")
endif()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/share/doc")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE.txt")
