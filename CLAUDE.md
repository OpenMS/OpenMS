# OpenMS Development

- Build folder: `OpenMS-build/`
- Build: `cmake --build OpenMS-build -j$(nproc)`
- Test: `ctest --test-dir OpenMS-build`

## Build Gotchas

- **Boost on macOS**: Use `-DBOOST_USE_STATIC_LIBS=OFF` (upstream Boost issue with transitive deps)
- **CMAKE_PREFIX_PATH**: Use `;` separator with `-D` flag, `:` with env var on Unix
- **delocate-wheel `-L`**: Destination subdir inside wheel, NOT library search path. Use `delocate-wheel --require-archs {delocate_archs} -w {dest_dir} -v {wheel}`

## Guidelines

- When CI tests fail, always investigate the root cause first before patching reference files or test expectations. Never bulk-replace expected values or regenerate test outputs as a first approach.
- Never dismiss a test failure as "acceptable behavior" or "pre-existing/unrelated" without explicit user confirmation. Always investigate fully before concluding a failure is not actionable.
- When asked to fix something, prefer using existing APIs/utilities in the codebase over heuristic approaches. Search for existing helper functions, conversion utilities, or detection APIs before implementing custom logic.
- Before implementing non-trivial changes, first explore the codebase for existing patterns and utilities, then explain the planned approach before writing code.
- For multi-file changes or unfamiliar areas, use plan mode or explore the relevant code first. Do not dive into implementation without understanding the existing architecture.

## pyOpenMS

See `src/pyOpenMS/CLAUDE.md` for Python binding development using nanobind.

Quick commands:
- Build: `cmake --build OpenMS-build --target pyopenms -j$(nproc)`
- Test: `PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/ -v`
