Custom Compilation of OpenMS
===========================

To compile with self built compilers and non default standard libraries, follow listed steps.

To choose any specific compiler, instead of the system default, add the whole path to these options for the cmake call:

```{tab} GCC
cmake -DCMAKE_C_COMPILER=/path/to/c-compiler/binary/gcc -DCMAKE_CXX_COMPILER=/path/to/c++-compiler/binary/g++
```

```{tab} Clang
cmake -DCMAKE_C_COMPILER=/path/to/c-compiler/binary/clang -DCMAKE_CXX_COMPILER=/path/to/c++-compiler/binary/clang++
```

To compile OpenMS with clang and a specific GCC stdlib, instead of the system default one:

Use this cmake option to specify an additional compiling option for clang:

```bash
cmake -DMY_CXX_FLAGS="--gcc-toolchain=/path/to/gcc"
```

with the path to the top gcc directory (containing the directory lib64) to the cmake call.

```{warning}
This combination does not work for all versions of clang and gcc.
- Clang 9.0.0 and GCC 4.8.5 stdlib does not work!
- Clang 9.0.0 and GCC 9.2.0 stdlib does not work!
- Clang 9.0.0 and GCC 8.3.0 stdlib compiles, but some tests fail.
- Clang 6.0.0 and GCC 7.4.0 stdlib (Ubuntu 18.04 default versions) works
```


