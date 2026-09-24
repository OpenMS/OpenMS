set(VCPKG_TARGET_ARCHITECTURE arm64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE dynamic)
set(VCPKG_BUILD_TYPE release)

set(VCPKG_CMAKE_SYSTEM_NAME Linux)

set(VCPKG_FIXUP_ELF_RPATH ON)

# vcpkg's openblas port otherwise lets OpenBLAS choose its kernels for the CPU of the
# build machine, and the binary cache passes that build on to later builds and to the
# packages; a build on a Neoverse N2 runner uses SVE, which other ARM CPUs lack. Build
# for a fixed target instead: ARMV8, which runs on every 64-bit ARM CPU.
if(PORT STREQUAL "openblas")
    list(APPEND VCPKG_CMAKE_CONFIGURE_OPTIONS "-DTARGET=ARMV8")
endif()
