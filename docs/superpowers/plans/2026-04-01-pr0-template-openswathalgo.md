# PR 0: Template OpenSwathAlgo for Float — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Template the openswathalgo scoring and stats helper functions to accept both `float` and `double`, and add a `vector<float>` overload to `IFeature::getIntensity()`.

**Architecture:** Template 11 Scoring functions and 3 StatsHelpers functions on scalar type `T`, using `Eigen::Matrix<T, Dynamic, 1>` internally. Add explicit template instantiations for `float` and `double` in `.cpp` files. Add `IFeature::getIntensity(vector<float>&)` as non-pure-virtual with default double-to-float conversion. Override in `FeatureOpenMS`.

**Tech Stack:** C++17, Eigen3, nanobind (no binding changes needed), OpenMS test framework (START_TEST/START_SECTION macros)

**Spec:** `docs/superpowers/specs/2026-04-01-float-intensity-propagation-design.md` (PR 0 section)

**Build:** `cmake --build OpenMS-build --target OpenSwathAlgo -j$(nproc)` for the library, `cmake --build OpenMS-build --target Scoring_test DiaHelpers_test -j$(nproc)` for tests

**Test:** `ctest --test-dir OpenMS-build -R "Scoring_test|DiaHelpers_test|StatsHelpers_test" -V`

---

### Task 1: Template `detail::normalize_sum_impl` and `normalize_sum`

**Files:**
- Modify: `src/openswathalgo/source/OPENSWATHALGO/ALGO/Scoring.cpp:19-38`
- Modify: `src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h:147,153`

- [ ] **Step 1: Template the internal helper in Scoring.cpp**

Change the `detail::normalize_sum_impl` from `double*` to a template, and template the `normalize_sum(vector)` overload. The deprecated `normalize_sum(double[], unsigned int)` stays as-is.

In `Scoring.cpp`, replace lines 19-38:

```cpp
namespace detail
{
  template <typename T>
  static void normalize_sum_impl(T* x, unsigned int n)
  {
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> v(x, n);
    T s = v.sum();
    if (s != T(0))
    {
      v /= s;
    }
  }
} // namespace detail

// deprecated C-array overload stays double
void normalize_sum(double x[], unsigned int n)
{
  detail::normalize_sum_impl<double>(x, n);
}

template <typename T>
void normalize_sum(std::vector<T>& x)
{
  detail::normalize_sum_impl<T>(x.data(), static_cast<unsigned int>(x.size()));
}

// explicit instantiations
template OPENSWATHALGO_DLLAPI void normalize_sum<double>(std::vector<double>&);
template OPENSWATHALGO_DLLAPI void normalize_sum<float>(std::vector<float>&);
```

In `Scoring.h`, replace the vector overload declaration at line 153:

```cpp
template <typename T>
OPENSWATHALGO_DLLAPI void normalize_sum(std::vector<T>& x);
```

Keep the deprecated C-array declaration at line 147 unchanged.

- [ ] **Step 2: Build and run existing tests**

```bash
cmake --build OpenMS-build --target OpenSwathAlgo Scoring_test -j$(nproc)
ctest --test-dir OpenMS-build -R Scoring_test -V
```

Expected: All existing tests pass (double path unchanged).

- [ ] **Step 3: Commit**

```bash
git add src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h \
        src/openswathalgo/source/OPENSWATHALGO/ALGO/Scoring.cpp
git commit -m "feat(openswathalgo): template normalize_sum for float/double"
```

---

### Task 2: Template `NormalizedManhattanDist`, `RootMeanSquareDeviation`, `SpectralAngle` (vector overloads)

**Files:**
- Modify: `src/openswathalgo/source/OPENSWATHALGO/ALGO/Scoring.cpp:51-109`
- Modify: `src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h:82,100,118`

- [ ] **Step 1: Template the three vector-overload functions**

In `Scoring.cpp`, replace `NormalizedManhattanDist(vector<double>&, vector<double>&)` (lines 51-59):

```cpp
template <typename T>
double NormalizedManhattanDist(std::vector<T>& x, std::vector<T>& y)
{
  OPENMS_PRECONDITION(x.size() == y.size(), "Vectors must be same size");
  // local copies for normalization
  std::vector<T> a(x), b(y);
  normalize_sum(a);
  normalize_sum(b);
  using VecT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  Eigen::Map<VecT> va(a.data(), a.size());
  Eigen::Map<VecT> vb(b.data(), b.size());
  return static_cast<double>((va - vb).cwiseAbs().sum());
}

template OPENSWATHALGO_DLLAPI double NormalizedManhattanDist<double>(std::vector<double>&, std::vector<double>&);
template OPENSWATHALGO_DLLAPI double NormalizedManhattanDist<float>(std::vector<float>&, std::vector<float>&);
```

Replace `RootMeanSquareDeviation(const vector<double>&, const vector<double>&)` (lines 70-76):

```cpp
template <typename T>
double RootMeanSquareDeviation(const std::vector<T>& x, const std::vector<T>& y)
{
  OPENMS_PRECONDITION(x.size() == y.size(), "Vectors must be same size");
  using VecT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  Eigen::Map<const VecT> va(x.data(), x.size());
  Eigen::Map<const VecT> vb(y.data(), y.size());
  return std::sqrt(static_cast<double>((va - vb).squaredNorm()) / x.size());
}

template OPENSWATHALGO_DLLAPI double RootMeanSquareDeviation<double>(const std::vector<double>&, const std::vector<double>&);
template OPENSWATHALGO_DLLAPI double RootMeanSquareDeviation<float>(const std::vector<float>&, const std::vector<float>&);
```

Replace `SpectralAngle(const vector<double>&, const vector<double>&)` (lines 98-109):

```cpp
template <typename T>
double SpectralAngle(const std::vector<T>& x, const std::vector<T>& y)
{
  OPENMS_PRECONDITION(x.size() == y.size(), "Vectors must be same size");
  using VecT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  Eigen::Map<const VecT> va(x.data(), x.size());
  Eigen::Map<const VecT> vb(y.data(), y.size());
  double x_len = static_cast<double>(va.norm());
  double y_len = static_cast<double>(vb.norm());
  double denominator = x_len * y_len;
  if (denominator < 1e-100) return 0.0;
  double theta = static_cast<double>(va.dot(vb)) / denominator;
  return std::acos(std::max(-1.0, std::min(1.0, theta)));
}

template OPENSWATHALGO_DLLAPI double SpectralAngle<double>(const std::vector<double>&, const std::vector<double>&);
template OPENSWATHALGO_DLLAPI double SpectralAngle<float>(const std::vector<float>&, const std::vector<float>&);
```

In `Scoring.h`, update the three vector-overload declarations (keep deprecated C-array versions unchanged):

```cpp
// Line 82:
template <typename T>
OPENSWATHALGO_DLLAPI double NormalizedManhattanDist(std::vector<T>& x, std::vector<T>& y);

// Line 100:
template <typename T>
OPENSWATHALGO_DLLAPI double RootMeanSquareDeviation(const std::vector<T>& x, const std::vector<T>& y);

// Line 118:
template <typename T>
OPENSWATHALGO_DLLAPI double SpectralAngle(const std::vector<T>& x, const std::vector<T>& y);
```

- [ ] **Step 2: Build and run tests**

```bash
cmake --build OpenMS-build --target OpenSwathAlgo Scoring_test -j$(nproc)
ctest --test-dir OpenMS-build -R Scoring_test -V
```

Expected: All existing tests pass.

- [ ] **Step 3: Commit**

```bash
git add src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h \
        src/openswathalgo/source/OPENSWATHALGO/ALGO/Scoring.cpp
git commit -m "feat(openswathalgo): template NormalizedManhattanDist, RMSD, SpectralAngle for float/double"
```

---

### Task 3: Template `standardize_data` and cross-correlation functions

**Files:**
- Modify: `src/openswathalgo/source/OPENSWATHALGO/ALGO/Scoring.cpp:128-261`
- Modify: `src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h:123-143`

- [ ] **Step 1: Template `standardize_data`**

Replace lines 128-148 in Scoring.cpp:

```cpp
template <typename T>
void standardize_data(std::vector<T>& data)
{
  if (data.empty()) return;
  using VecT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  Eigen::Map<VecT> v(data.data(), data.size());
  T mean = v.mean();
  T stdev = std::sqrt((v.array() - mean).square().mean());
  if (mean == T(0) && stdev == T(0))
  {
    return;
  }
  if (stdev == T(0))
  {
    v.setZero();
    return;
  }
  v = (v.array() - mean) / stdev;
}

template OPENSWATHALGO_DLLAPI void standardize_data<double>(std::vector<double>&);
template OPENSWATHALGO_DLLAPI void standardize_data<float>(std::vector<float>&);
```

- [ ] **Step 2: Template `calculateCrossCorrelation`**

Replace lines 173-200 in Scoring.cpp. Key changes: `const double* __restrict` becomes `const T* __restrict`, `double sxy` stays double (feeds into `XCorrEntry.second` which is double), Eigen maps use `Matrix<T, Dynamic, 1>`:

```cpp
template <typename T>
XCorrArrayType calculateCrossCorrelation(const std::vector<T>& data1,
    const std::vector<T>& data2, const int maxdelay, const int lag)
{
  XCorrArrayType result;
  const int datasize = static_cast<int>(data1.size());
  if (datasize == 0) return result;
  OPENMS_PRECONDITION(datasize == static_cast<int>(data2.size()), "data1 and data2 must have same size");

  const T* __restrict d1 = data1.data();
  const T* __restrict d2 = data2.data();

  for (int delay = -maxdelay; delay <= maxdelay; delay += lag)
  {
    double sxy = 0.0;
    int start = std::max(0, delay);
    int end = std::min(datasize, datasize + delay);
    int len = end - start;
    if (len > 0)
    {
      using VecT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
      Eigen::Map<const VecT> sub1(d1 + start, len);
      Eigen::Map<const VecT> sub2(d2 + start - delay, len);
      sxy = static_cast<double>(sub1.dot(sub2));
    }
    result.data.emplace_back(delay, sxy);
  }
  return result;
}

template OPENSWATHALGO_DLLAPI XCorrArrayType calculateCrossCorrelation<double>(const std::vector<double>&, const std::vector<double>&, const int, const int);
template OPENSWATHALGO_DLLAPI XCorrArrayType calculateCrossCorrelation<float>(const std::vector<float>&, const std::vector<float>&, const int, const int);
```

- [ ] **Step 3: Template `normalizedCrossCorrelationPost` and `normalizedCrossCorrelation`**

Replace lines 150-171:

```cpp
template <typename T>
XCorrArrayType normalizedCrossCorrelation(std::vector<T>& data1,
    std::vector<T>& data2, const int maxdelay, const int lag)
{
  standardize_data(data1);
  standardize_data(data2);
  return normalizedCrossCorrelationPost(data1, data2, maxdelay, lag);
}

template <typename T>
XCorrArrayType normalizedCrossCorrelationPost(std::vector<T>& normalized_data1,
    std::vector<T>& normalized_data2, const int maxdelay, const int lag)
{
  auto result = calculateCrossCorrelation(normalized_data1, normalized_data2, maxdelay, lag);
  for (auto& it : result.data)
  {
    it.second /= static_cast<double>(normalized_data1.size());
  }
  return result;
}

// explicit instantiations for all four functions
template OPENSWATHALGO_DLLAPI XCorrArrayType normalizedCrossCorrelation<double>(std::vector<double>&, std::vector<double>&, const int, const int);
template OPENSWATHALGO_DLLAPI XCorrArrayType normalizedCrossCorrelation<float>(std::vector<float>&, std::vector<float>&, const int, const int);
template OPENSWATHALGO_DLLAPI XCorrArrayType normalizedCrossCorrelationPost<double>(std::vector<double>&, std::vector<double>&, const int, const int);
template OPENSWATHALGO_DLLAPI XCorrArrayType normalizedCrossCorrelationPost<float>(std::vector<float>&, std::vector<float>&, const int, const int);
```

- [ ] **Step 4: Template `calcxcorr_legacy_mquest_`**

Replace lines 202-261. Change all `Eigen::Map<Eigen::VectorXd>` to `Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>`, all `double` internal variables (`mean1`, `mean2`, `sqsum1`, `sqsum2`, `sxy`, `denominator`) to `T`, and `double sxy = 0.0` to `T sxy = T(0)`. The return `XCorrEntry` second element stays double via `static_cast<double>(sxy / denominator)`.

```cpp
template <typename T>
XCorrArrayType calcxcorr_legacy_mquest_(std::vector<T>& data1,
    std::vector<T>& data2, bool normalize)
{
  // [full implementation with T replacing double for internal vars and Eigen maps]
  // ... (follow existing structure, change types)
}

template OPENSWATHALGO_DLLAPI XCorrArrayType calcxcorr_legacy_mquest_<double>(std::vector<double>&, std::vector<double>&, bool);
template OPENSWATHALGO_DLLAPI XCorrArrayType calcxcorr_legacy_mquest_<float>(std::vector<float>&, std::vector<float>&, bool);
```

- [ ] **Step 5: Update Scoring.h declarations**

Update lines 123-143 to templated declarations:

```cpp
template <typename T>
OPENSWATHALGO_DLLAPI XCorrArrayType calcxcorr_legacy_mquest_(std::vector<T>& data1,
    std::vector<T>& data2, bool normalize);

template <typename T>
OPENSWATHALGO_DLLAPI XCorrArrayType normalizedCrossCorrelation(std::vector<T>& data1,
    std::vector<T>& data2, const int maxdelay, const int lag);

template <typename T>
OPENSWATHALGO_DLLAPI XCorrArrayType normalizedCrossCorrelationPost(std::vector<T>& normalized_data1,
    std::vector<T>& normalized_data2, const int maxdelay, const int lag);

template <typename T>
OPENSWATHALGO_DLLAPI XCorrArrayType calculateCrossCorrelation(const std::vector<T>& data1,
    const std::vector<T>& data2, const int maxdelay, const int lag);

// xcorrArrayGetMaxPeak stays non-templated (takes XCorrArrayType, not vector<T>)

template <typename T>
OPENSWATHALGO_DLLAPI void standardize_data(std::vector<T>& data);
```

- [ ] **Step 6: Build and run tests**

```bash
cmake --build OpenMS-build --target OpenSwathAlgo Scoring_test -j$(nproc)
ctest --test-dir OpenMS-build -R Scoring_test -V
```

Expected: All existing tests pass (double instantiation used).

- [ ] **Step 7: Commit**

```bash
git add src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h \
        src/openswathalgo/source/OPENSWATHALGO/ALGO/Scoring.cpp
git commit -m "feat(openswathalgo): template cross-correlation and standardize_data for float/double"
```

---

### Task 4: Template `computeAndAppendRank` and `computeRankVector`

**Files:**
- Modify: `src/openswathalgo/source/OPENSWATHALGO/ALGO/Scoring.cpp:263-295`
- Modify: `src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h:156,159`

- [ ] **Step 1: Template the rank functions**

In Scoring.cpp, replace lines 263-295:

```cpp
template <typename T>
void computeAndAppendRank(const std::vector<T>& v_temp, std::vector<unsigned int>& ranks)
{
  std::vector<unsigned int> indices(v_temp.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
    [&v_temp](unsigned int a, unsigned int b) { return v_temp[a] < v_temp[b]; });

  unsigned int rank = 0;
  T x = T(0);
  for (unsigned int i = 0; i < indices.size(); ++i)
  {
    if (i == 0 || v_temp[indices[i]] != x)
    {
      rank = i + 1;
    }
    x = v_temp[indices[i]];
    ranks.push_back(rank);
  }
}

template <typename T>
std::vector<unsigned int> computeRankVector(const std::vector<std::vector<T>>& intensity,
    std::vector<std::vector<unsigned int>>& ranks)
{
  ranks.clear();
  ranks.resize(intensity.size());
  for (Size i = 0; i < intensity.size(); ++i)
  {
    computeAndAppendRank(intensity[i], ranks[i]);
  }
  // [rest of existing logic for combined rank computation stays unchanged]
  std::vector<unsigned int> result;
  // ... existing code
  return result;
}

template OPENSWATHALGO_DLLAPI void computeAndAppendRank<double>(const std::vector<double>&, std::vector<unsigned int>&);
template OPENSWATHALGO_DLLAPI void computeAndAppendRank<float>(const std::vector<float>&, std::vector<unsigned int>&);
template OPENSWATHALGO_DLLAPI std::vector<unsigned int> computeRankVector<double>(const std::vector<std::vector<double>>&, std::vector<std::vector<unsigned int>>&);
template OPENSWATHALGO_DLLAPI std::vector<unsigned int> computeRankVector<float>(const std::vector<std::vector<float>>&, std::vector<std::vector<unsigned int>>&);
```

Update Scoring.h declarations:

```cpp
template <typename T>
OPENSWATHALGO_DLLAPI void computeAndAppendRank(const std::vector<T>& v_temp, std::vector<unsigned int>& ranks);

template <typename T>
OPENSWATHALGO_DLLAPI std::vector<unsigned int> computeRankVector(const std::vector<std::vector<T>>& intensity,
    std::vector<std::vector<unsigned int>>& ranks);
```

- [ ] **Step 2: Build and run tests**

```bash
cmake --build OpenMS-build --target OpenSwathAlgo Scoring_test -j$(nproc)
ctest --test-dir OpenMS-build -R Scoring_test -V
```

- [ ] **Step 3: Commit**

```bash
git add src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h \
        src/openswathalgo/source/OPENSWATHALGO/ALGO/Scoring.cpp
git commit -m "feat(openswathalgo): template computeAndAppendRank/computeRankVector for float/double"
```

---

### Task 5: Template StatsHelpers (`normalize`, `dotprodScoring`, `manhattanScoring`) and add float instantiations for `norm`/`manhattanDist`

**Files:**
- Modify: `src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h:25,71,97`
- Modify: `src/openswathalgo/source/OPENSWATHALGO/ALGO/StatsHelpers.cpp:20-110`

- [ ] **Step 1: Template `normalize`**

In StatsHelpers.cpp, replace lines 20-36:

```cpp
template <typename T>
void normalize(const std::vector<T>& intensities, double normalization_factor,
               std::vector<T>& normalized_intensities)
{
  normalized_intensities.resize(intensities.size());
  if (normalization_factor != 0.0)
  {
    T factor = static_cast<T>(1.0 / normalization_factor);
    std::transform(intensities.begin(), intensities.end(),
                   normalized_intensities.begin(),
                   [factor](T val) { return val * factor; });
  }
  else
  {
    std::fill(normalized_intensities.begin(), normalized_intensities.end(), T(0));
  }
}

template OPENSWATHALGO_DLLAPI void normalize<double>(const std::vector<double>&, double, std::vector<double>&);
template OPENSWATHALGO_DLLAPI void normalize<float>(const std::vector<float>&, double, std::vector<float>&);
```

In StatsHelpers.h, update line 25:

```cpp
template <typename T>
OPENSWATHALGO_DLLAPI void normalize(const std::vector<T>& intensities,
    double normalization_factor, std::vector<T>& normalized_intensities);
```

- [ ] **Step 2: Template `dotprodScoring` and `manhattanScoring`**

In StatsHelpers.cpp, replace lines 39-72:

```cpp
template <typename T>
double dotprodScoring(std::vector<T> intExp, std::vector<T> theorint)
{
  using VecT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  Eigen::Map<VecT> mapExp(intExp.data(), intExp.size());
  Eigen::Map<VecT> mapTheor(theorint.data(), theorint.size());
  mapExp = mapExp.cwiseSqrt();
  mapTheor = mapTheor.cwiseSqrt();

  double expNorm = norm(intExp.begin(), intExp.end());
  normalize(intExp, expNorm, intExp);
  double theoNorm = norm(theorint.begin(), theorint.end());
  normalize(theorint, theoNorm, theorint);

  return dotProd(intExp.begin(), intExp.end(), theorint.begin());
}

template <typename T>
double manhattanScoring(std::vector<T> intExp, std::vector<T> theorint)
{
  using VecT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  Eigen::Map<VecT> mapExp(intExp.data(), intExp.size());
  Eigen::Map<VecT> mapTheor(theorint.data(), theorint.size());
  mapExp = mapExp.cwiseSqrt();
  mapTheor = mapTheor.cwiseSqrt();

  double expSum = static_cast<double>(mapExp.sum());
  double theoSum = static_cast<double>(mapTheor.sum());
  normalize(intExp, expSum, intExp);
  normalize(theorint, theoSum, theorint);

  return manhattanDist(intExp.begin(), intExp.end(), theorint.begin());
}

template OPENSWATHALGO_DLLAPI double dotprodScoring<double>(std::vector<double>, std::vector<double>);
template OPENSWATHALGO_DLLAPI double dotprodScoring<float>(std::vector<float>, std::vector<float>);
template OPENSWATHALGO_DLLAPI double manhattanScoring<double>(std::vector<double>, std::vector<double>);
template OPENSWATHALGO_DLLAPI double manhattanScoring<float>(std::vector<float>, std::vector<float>);
```

Update StatsHelpers.h:

```cpp
template <typename T>
OPENSWATHALGO_DLLAPI double dotprodScoring(std::vector<T> intExp, std::vector<T> theorint);

template <typename T>
OPENSWATHALGO_DLLAPI double manhattanScoring(std::vector<T> intExp, std::vector<T> theorint);
```

- [ ] **Step 3: Add float iterator instantiations for `norm` and `manhattanDist`**

In StatsHelpers.cpp, add after the existing double instantiations:

After the existing `norm` instantiations (around line 90), add:

```cpp
template OPENSWATHALGO_DLLAPI double norm<std::vector<float>::const_iterator>(std::vector<float>::const_iterator, std::vector<float>::const_iterator);
template OPENSWATHALGO_DLLAPI double norm<std::vector<float>::iterator>(std::vector<float>::iterator, std::vector<float>::iterator);
```

After the existing `manhattanDist` instantiations (around line 110), add:

```cpp
template OPENSWATHALGO_DLLAPI double manhattanDist<std::vector<float>::const_iterator, std::vector<float>::const_iterator>(std::vector<float>::const_iterator, std::vector<float>::const_iterator, std::vector<float>::const_iterator);
```

Also add the corresponding `extern template` declarations in StatsHelpers.h alongside the existing double ones.

- [ ] **Step 4: Build and run tests**

```bash
cmake --build OpenMS-build --target OpenSwathAlgo Scoring_test DiaHelpers_test -j$(nproc)
ctest --test-dir OpenMS-build -R "Scoring_test|DiaHelpers_test|StatsHelpers_test" -V
```

Expected: All existing tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h \
        src/openswathalgo/source/OPENSWATHALGO/ALGO/StatsHelpers.cpp
git commit -m "feat(openswathalgo): template normalize, dotprodScoring, manhattanScoring for float/double"
```

---

### Task 6: Add float test cases to `Scoring_test.cpp`

**Files:**
- Modify: `src/tests/class_tests/openswathalgo/Scoring_test.cpp`

- [ ] **Step 1: Add float cross-correlation test section**

Add after the existing `normalizedCrossCorrelation` test section (around line 300):

```cpp
START_SECTION(normalizedCrossCorrelation_float)
{
  // Same data as double test, but using vector<float>
  std::vector<float> data1_f = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
  std::vector<float> data2_f = {5.0f, 4.0f, 3.0f, 2.0f, 1.0f};

  auto result = OpenSwath::Scoring::normalizedCrossCorrelation(data1_f, data2_f, 2, 1);
  // Verify structure
  TEST_EQUAL(result.data.size(), 5)

  // Compare against double result with tolerance
  std::vector<double> data1_d = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> data2_d = {5.0, 4.0, 3.0, 2.0, 1.0};
  auto result_d = OpenSwath::Scoring::normalizedCrossCorrelation(data1_d, data2_d, 2, 1);

  for (Size i = 0; i < result.data.size(); ++i)
  {
    TEST_EQUAL(result.data[i].first, result_d.data[i].first)
    TEST_REAL_SIMILAR(result.data[i].second, result_d.data[i].second)
  }
}
END_SECTION
```

- [ ] **Step 2: Add float rank computation test section**

```cpp
START_SECTION(computeRankVector_float)
{
  std::vector<std::vector<float>> intensity_f = {{3.0f, 1.0f, 2.0f}, {6.0f, 4.0f, 5.0f}};
  std::vector<std::vector<unsigned int>> ranks;
  auto combined = OpenSwath::Scoring::computeRankVector(intensity_f, ranks);
  TEST_EQUAL(ranks.size(), 2)
  TEST_EQUAL(ranks[0].size(), 3)
  // Rank order should match double version
  TEST_EQUAL(ranks[0][0], 3) // 3.0 is rank 3
  TEST_EQUAL(ranks[0][1], 1) // 1.0 is rank 1
  TEST_EQUAL(ranks[0][2], 2) // 2.0 is rank 2
}
END_SECTION
```

- [ ] **Step 3: Add float standardize_data test**

```cpp
START_SECTION(standardize_data_float)
{
  std::vector<float> data_f = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
  std::vector<double> data_d = {1.0, 2.0, 3.0, 4.0, 5.0};
  OpenSwath::Scoring::standardize_data(data_f);
  OpenSwath::Scoring::standardize_data(data_d);
  for (Size i = 0; i < data_f.size(); ++i)
  {
    TEST_REAL_SIMILAR(static_cast<double>(data_f[i]), data_d[i])
  }
}
END_SECTION
```

- [ ] **Step 4: Build and run**

```bash
cmake --build OpenMS-build --target Scoring_test -j$(nproc)
ctest --test-dir OpenMS-build -R Scoring_test -V
```

Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/tests/class_tests/openswathalgo/Scoring_test.cpp
git commit -m "test(openswathalgo): add float instantiation tests for Scoring functions"
```

---

### Task 7: Add float test cases to `DiaHelpers_test.cpp`

**Files:**
- Modify: `src/tests/class_tests/openswathalgo/DiaHelpers_test.cpp`

- [ ] **Step 1: Add float dotprodScoring and manhattanScoring tests**

Add a new section after the existing test:

```cpp
START_SECTION(testDotProdScore_float)
{
  float arr1[] = {100.0f, 200.0f, 1.0f, 50.0f, 300.0f, 2.0f};
  float arr2[] = {200.0f, 100.0f, 2.0f, 300.0f, 50.0f, 1.0f};
  std::vector<float> intExp(arr1, arr1 + 6);
  std::vector<float> theorint(arr2, arr2 + 6);

  double dotprod = OpenSwath::dotprodScoring(intExp, theorint);
  TEST_REAL_SIMILAR(dotprod, 0.8604286) // same expected value

  double manhattan = OpenSwath::manhattanScoring(intExp, theorint);
  TEST_REAL_SIMILAR(manhattan, 0.4950837)
}
END_SECTION
```

- [ ] **Step 2: Build and run**

```bash
cmake --build OpenMS-build --target DiaHelpers_test -j$(nproc)
ctest --test-dir OpenMS-build -R DiaHelpers_test -V
```

- [ ] **Step 3: Commit**

```bash
git add src/tests/class_tests/openswathalgo/DiaHelpers_test.cpp
git commit -m "test(openswathalgo): add float instantiation tests for dotprodScoring/manhattanScoring"
```

---

### Task 8: Add `IFeature::getIntensity(vector<float>&)` overload

**Files:**
- Modify: `src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h:20-28`
- Modify: `src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/MockObjects.h`
- Modify: `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h`
- Modify: `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.cpp`

- [ ] **Step 1: Add non-pure-virtual float overload to IFeature**

In `ITransition.h`, add after line 25 (`virtual void getIntensity(std::vector<double>& intens) const = 0;`):

```cpp
    virtual void getIntensity(std::vector<float>& intens) const
    {
      // Default implementation: delegate to double version and narrow
      std::vector<double> tmp;
      getIntensity(tmp);
      intens.assign(tmp.begin(), tmp.end());
    }
```

- [ ] **Step 2: Add float overload to MockFeature**

In `MockObjects.h`, add a member variable and override inside class `MockFeature`:

```cpp
    void getIntensity(std::vector<float>& intens) const override
    {
      intens.assign(m_intensity_vec.begin(), m_intensity_vec.end());
    }
```

This converts the existing `std::vector<double> m_intensity_vec` to float on demand.

- [ ] **Step 3: Add float overload to FeatureOpenMS**

In `MRMFeatureAccessOpenMS.h`, add declaration in class `FeatureOpenMS`:

```cpp
    void getIntensity(std::vector<float>& intens) const override;
```

In `MRMFeatureAccessOpenMS.cpp`, add implementation after the existing `getIntensity(vector<double>&)`:

```cpp
void FeatureOpenMS::getIntensity(std::vector<float>& intens) const
{
  OPENMS_PRECONDITION(feature_->getConvexHulls().size() == 1, "Must have exactly one convex hull");
  intens.clear();
  const auto& points = feature_->getConvexHulls()[0].getHullPoints();
  intens.reserve(points.size());
  for (const auto& pt : points)
  {
    intens.push_back(static_cast<float>(pt.getY()));
  }
}
```

- [ ] **Step 4: Build the full project and run tests**

```bash
cmake --build OpenMS-build --target OpenSwathAlgo pyopenms -j$(nproc) 2>&1 | head -50
ctest --test-dir OpenMS-build -R "Scoring_test|DiaHelpers_test|OpenSwathMRMFeatureAccessOpenMS_test" -V
```

Expected: All tests pass. No compilation errors.

- [ ] **Step 5: Commit**

```bash
git add src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h \
        src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/MockObjects.h \
        src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h \
        src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.cpp
git commit -m "feat(openswathalgo): add IFeature::getIntensity(vector<float>&) overload"
```

---

### Task 9: Final integration test — build everything and run full test suite

**Files:** None (verification only)

- [ ] **Step 1: Full build**

```bash
cmake --build OpenMS-build -j$(nproc)
```

Expected: Clean build with no warnings about the changed files.

- [ ] **Step 2: Run all OpenSWATH-related tests**

```bash
ctest --test-dir OpenMS-build -R "Scoring_test|DiaHelpers_test|StatsHelpers_test|MRMScoring_test|OpenSwathScoring_test|MRMFeatureFinderScoring_test|OpenSwathMRMFeatureAccessOpenMS" -V
```

Expected: All tests pass. The double instantiation is used throughout — no behavioral change yet.

- [ ] **Step 3: Verify float instantiations link correctly**

Check that the float template instantiations are exported from the shared library:

```bash
nm -C OpenMS-build/lib/libOpenSwathAlgo.so | grep "normalize_sum.*float" | head -3
nm -C OpenMS-build/lib/libOpenSwathAlgo.so | grep "standardize_data.*float" | head -3
nm -C OpenMS-build/lib/libOpenSwathAlgo.so | grep "normalizedCrossCorrelation.*float" | head -3
```

Expected: Symbol entries for all float instantiations.

- [ ] **Step 4: Final commit (if any fixups needed)**

If any test failures or build warnings were found and fixed:

```bash
git add -u
git commit -m "fix(openswathalgo): address build/test issues in float templating"
```
