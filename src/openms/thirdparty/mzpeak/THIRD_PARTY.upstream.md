# Third-party code

This project vendors a small amount of third-party source. Each retains its
upstream copyright and license header.

## ms-numpress (MSNumpress)

- **Files:** `src/util/vendor/MSNumpress.cpp`, `include/mzpeak/util/vendor/MSNumpress.hpp`
- **Upstream:** https://github.com/ms-numpress/ms-numpress
- **Copyright:** 2013 Johan Teleman
- **License:** Apache License 2.0 (see the header block in each file)
- **Why vendored:** numpress is the codec behind the chunked mzPeak transforms
  `MS:1002312` (linear, m/z) and `MS:1002314` (SLOF, intensity). None of the
  build dependencies (Arrow/Parquet, Boost, libzip) provide it, and it is not
  packaged for meson/nix, so the two-file Apache-2.0 implementation is copied in.
- **Local modifications:** none of substance — only the `#include` path of the
  header was changed to `"mzpeak/util/vendor/MSNumpress.hpp"` so it resolves
  inside this include tree. The file is intentionally excluded from
  `clang-format` to stay diffable against upstream. A thin wrapper
  (`src/util/numpress.cpp`, `include/mzpeak/util/numpress.h`) adapts the API to
  this project's namespace and types.
