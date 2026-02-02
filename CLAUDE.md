# OpenMS Development

- Build folder: `OpenMS-build/`
- Build: `cmake --build OpenMS-build -j$(nproc)`
- Test: `ctest --test-dir OpenMS-build`

## Build Gotchas

- **Boost on macOS**: Use `-DBOOST_USE_STATIC_LIBS=OFF` (upstream Boost issue with transitive deps)
- **CMAKE_PREFIX_PATH**: Use `;` separator with `-D` flag, `:` with env var on Unix
- **delocate-wheel `-L`**: Destination subdir inside wheel, NOT library search path. Use `delocate-wheel --require-archs {delocate_archs} -w {dest_dir} -v {wheel}`

## pyOpenMS

See `src/pyOpenMS/CLAUDE.md` for Python binding development using cython.
See `src/pyOpenMS2/CLAUDE.md` for experimental Python binding development using nanobind.

Quick commands:
- Build: `cmake --build OpenMS-build --target pyopenms -j$(nproc)`
- Test: `PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/ -v`
- Force rebuild after addon changes: `rm OpenMS-build/pyOpenMS/.cpp_extension_generated`
