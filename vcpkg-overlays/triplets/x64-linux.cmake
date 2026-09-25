set(VCPKG_TARGET_ARCHITECTURE x64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE static)

set(VCPKG_CMAKE_SYSTEM_NAME Linux)

# vcpkg's openblas port otherwise lets OpenBLAS choose its kernels for the CPU of the
# build machine, and the binary cache passes that build on to later builds and to the
# packages. Build for a fixed target instead: CORE2 (SSSE3), the x86_64 baseline OpenMS
# itself is compiled for (-mssse3), so the packages run on every CPU OpenMS supports.
# OpenMS calls BLAS only through COIN-OR's LAPACK, so faster kernels would gain it
# little. The fixed target also keeps clang debug builds working on AVX-512 machines,
# where OpenBLAS's AVX-512 SGEMM kernel does not compile without optimization.
#
# This is the default host triplet: vcpkg builds openblas for it as well, for the
# getarch tools the target triplet's openblas build runs, so it needs the same target.
if(PORT STREQUAL "openblas")
    list(APPEND VCPKG_CMAKE_CONFIGURE_OPTIONS "-DTARGET=CORE2")
endif()
