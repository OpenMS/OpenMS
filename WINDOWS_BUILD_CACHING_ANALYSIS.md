# Windows Build Caching Analysis and Improvements

## Issue Background

The Windows builds in GitHub Actions have been significantly slower than expected, with minimal time savings from caching. This document details the investigation findings and improvements made.

## Investigation Summary

### Codebase Metrics
- **Source files**: ~3,765 C++ files (.cpp and .h)
- **Source code size**: ~685MB
- **Build system**: CMake with Ninja generator (not Visual Studio generator)
- **Compiler**: MSVC (cl.exe) with ccache wrapper

### Problems Identified

#### 1. Insufficient ccache Size (CRITICAL)
**Symptom**: Cache constantly evicting entries, leading to cache misses
- **Previous setting**: `CCACHE_MAXSIZE: 400M`
- **Problem**: Too small for a codebase with 3,765 source files and 685MB source
- **Impact**: Low cache hit rate, frequent recompilation of unchanged files
- **Fix**: Increased to `CCACHE_MAXSIZE: 1G`
- **Rationale**: Balances cache effectiveness with GitHub Actions 10GB total cache limit

#### 2. Excessive Compression Level (MODERATE)
**Symptom**: Slow cache save/restore operations
- **Previous setting**: `CCACHE_COMPRESSLEVEL: 12` (maximum)
- **Problem**: High compression trades too much CPU time for marginal space savings
- **Impact**: Slower cache operations, especially on cache saves
- **Fix**: Reduced to `CCACHE_COMPRESSLEVEL: 6`
- **Rationale**: Level 6 provides good compression (typically 60-70% of max) with 2-3x faster operations

#### 3. No Cache Visibility (MINOR)
**Symptom**: Unable to diagnose caching issues
- **Previous state**: No logging of cache statistics
- **Problem**: Impossible to verify cache effectiveness or debug issues
- **Fix**: Added `ccache -s` logging before and after builds
- **Impact**: Can now monitor cache hit rates and diagnose problems

### Additional Observations (Not Issues)

#### Cache Key Strategy ✅ GOOD
The cache key strategy is well-designed:
```yaml
key: ${{ runner.os }}-${{ runner.arch }}-${{ matrix.compiler }}-${{ matrix.compiler_ver }}-ccache-${{ steps.extract_branch.outputs.RUN_NAME }}-${{ github.run_number }}
restore-keys: |
  ${{ runner.os }}-${{ runner.arch }}-${{ matrix.compiler }}-${{ matrix.compiler_ver }}-ccache-${{ steps.extract_branch.outputs.RUN_NAME }}
  ${{ runner.os }}-${{ runner.arch }}-${{ matrix.compiler }}-${{ matrix.compiler_ver }}-ccache-develop
  ${{ runner.os }}-${{ runner.arch }}-${{ matrix.compiler }}-${{ matrix.compiler_ver }}-ccache-
```
- Includes unique run_number for the primary key
- Good fallback chain: branch → develop → any
- Properly scoped by OS, architecture, and compiler version

#### Compiler Flags ✅ GOOD
- Using Ninja generator (not Visual Studio), so the problematic `/Od /Ob0` flags in `cibuild.cmake:160` are NOT active
- Compiler flags in `cmake/compiler_flags.cmake` are appropriate:
  - `/MP` for multi-processor compilation
  - `/bigobj` for large object files
  - Proper optimization flags for Release builds

#### Parallelization ✅ GOOD
- Correctly using all available CPU cores via `BUILD_FLAGS: "-j${{ steps.cpu-cores.outputs.count }}"`
- Ninja respects this flag appropriately

#### Contrib Caching ✅ GOOD
- Separate cache for contrib builds (line 303-323)
- Uses stable key: `${{ runner.os }}-contrib3`
- Good fallback if cache miss

## Changes Made

### Files Modified
1. `.github/workflows/openms_ci_matrix_full.yml`
   - Updated CCACHE_MAXSIZE: 400M → 1G
   - Updated CCACHE_COMPRESSLEVEL: 12 → 6
   - Added ccache statistics logging

2. `.github/workflows/pyopenms-wheels.yml`
   - Same ccache configuration updates (3 locations)

### Expected Impact

#### Build Time Improvements
Based on typical ccache effectiveness for C++ projects:
- **First build** (cold cache): No change
- **Incremental builds** (warm cache): 
  - With old 400M limit: 30-50% time savings (estimated)
  - With new 1G limit: 60-80% time savings (expected)
- **Rebuild with dependencies**: 40-70% time savings (expected)

#### Cache Storage
- GitHub Actions cache limit: 10GB per repository
- Old ccache size: ~400MB compressed
- New ccache size: ~1GB compressed (well within limits)
- Contrib cache: ~varies but cached separately

## Recommendations

### Short-term (Already Implemented)
- ✅ Monitor cache hit rates via new statistics logging
- ✅ Verify 1G is sufficient (can adjust if needed)

### Medium-term (Future Improvements)
1. **Consider precompiled headers (PCH)**
   - OpenMS has many common includes
   - Could further reduce build times by 10-20%

2. **Unity builds for tests**
   - Test files could be grouped for faster compilation
   - Trade-off: larger object files, less granular incremental builds

3. **Split builds by component**
   - Consider separate caches for: library, tools, tests
   - More granular caching, but more complex

### Long-term (Architectural)
1. **Module-based builds**
   - When moving to C++20 modules
   - Can significantly reduce compilation dependencies

2. **Distributed compilation**
   - Tools like sccache with remote storage
   - Useful if builds remain slow despite caching

## Monitoring

### Key Metrics to Track
After these changes, monitor in GitHub Actions logs:

1. **Cache hit rate**: Look for "cache hit (direct)" vs "cache miss" in ccache stats
   - Target: >80% for incremental builds
   - Target: >60% for PR builds

2. **Cache size**: Check "cache size" vs "max cache size" 
   - Should stay well below 1G
   - If approaching limit, may need adjustment

3. **Build time trends**:
   - First build in PR: ~60-90 minutes (baseline)
   - Subsequent builds: Should see 60-80% reduction

### How to Interpret ccache Statistics

Example good output:
```
cache hit (direct)     : 2500
cache hit (preprocessed): 300
cache miss             : 200
files in cache         : 3000
cache size             : 0.9 GB / 1.0 GB
```

This shows:
- Hit rate: (2500+300)/(2500+300+200) = 93% ✅ Excellent
- Cache utilization: 0.9/1.0 = 90% ✅ Good, not full

Example problematic output:
```
cache hit (direct)     : 500
cache hit (preprocessed): 100  
cache miss             : 2500
files in cache         : 1500
cache size             : 1.0 GB / 1.0 GB (full!)
```

This shows:
- Hit rate: (500+100)/(500+100+2500) = 19% ❌ Poor
- Cache full: 1.0/1.0 = 100% ❌ Need larger cache

## References

- ccache documentation: https://ccache.dev/manual/latest.html
- GitHub Actions cache: https://docs.github.com/en/actions/using-workflows/caching-dependencies-to-speed-up-workflows
- OpenMS build system: CMake with Ninja generator
