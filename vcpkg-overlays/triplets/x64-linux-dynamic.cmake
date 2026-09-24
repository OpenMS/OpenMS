set(VCPKG_TARGET_ARCHITECTURE x64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE dynamic)

set(VCPKG_CMAKE_SYSTEM_NAME Linux)

set(VCPKG_FIXUP_ELF_RPATH ON)

# vcpkg's openblas port otherwise lets OpenBLAS choose its kernels for the CPU of the
# build machine, and the binary cache passes that build on to later builds and to the
# packages. Build for a fixed target instead: HASWELL (AVX2). That also keeps clang
# debug builds working on AVX-512 machines, where OpenBLAS's AVX-512 SGEMM kernel
# does not compile without optimization.
if(PORT STREQUAL "openblas")
    list(APPEND VCPKG_CMAKE_CONFIGURE_OPTIONS "-DTARGET=HASWELL")
endif()
