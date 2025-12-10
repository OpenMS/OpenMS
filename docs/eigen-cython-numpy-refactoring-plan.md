# Refactoring Plan: Eigen-Cython-NumPy Interface

## Problem Statement

OpenMS currently exports Eigen symbols in its public API/ABI, primarily because:

1. **`Matrix<T>` inherits from `Eigen::Matrix<T>`** (`src/openms/include/OpenMS/DATASTRUCTURES/Matrix.h:31`)
2. **Public APIs expose Eigen types directly** (e.g., `TraceFitter`, `BinnedSpectrum`, `NonNegativeLeastSquaresSolver`)

This causes:
- ABI instability when Eigen is updated
- Larger symbol tables and binary sizes
- Tight coupling between OpenMS and Eigen versions
- Unnecessary complexity in Python bindings

## Goal

Hide Eigen from the public API while:
- Maintaining full Eigen functionality for internal C++ code
- Preserving efficient NumPy interoperability in Cython
- Minimizing breaking changes

---

## Current Architecture Analysis

### Classes Exposing Eigen in Public API

| Class | File | Exposed Types |
|-------|------|---------------|
| `Matrix<T>` | `DATASTRUCTURES/Matrix.h:31` | Inherits from `Eigen::Matrix` |
| `TraceFitter::GenericFunctor` | `FEATUREFINDER/TraceFitter.h:53-55` | `Eigen::VectorXd`, `Eigen::MatrixXd` |
| `BinnedSpectrum` | `KERNEL/BinnedSpectrum.h:87` | `Eigen::SparseVector<float>` |
| `NonNegativeLeastSquaresSolver` | `ML/NNLS/NonNegativeLeastSquaresSolver.h:52` | `Matrix<double>::EigenMatrixType` |
| `LevMarqFitter1D` | `FEATUREFINDER/LevMarqFitter1D.h` | `Eigen::VectorXd`, `Eigen::MatrixXd` |
| `EGHTraceFitter` | `FEATUREFINDER/EGHTraceFitter.h` | Via `TraceFitter` |
| `GaussTraceFitter` | `FEATUREFINDER/GaussTraceFitter.h` | Via `TraceFitter` |

### How Cython Currently Interfaces with NumPy

The Cython bindings (`src/pyOpenMS/addons/MatrixDouble.pyx`) only need:

```cython
cdef double* data = mat_.data()      # Raw pointer
cdef unsigned int rows = mat_.rows() # Dimensions
cdef unsigned int cols = mat_.cols()
cdef int inner = mat_.innerStride()  # Strides
cdef int outer = mat_.outerStride()
cdef bool row_major = mat_.rowMajor() # Layout
```

**Key insight: Cython does not need Eigen types—only raw buffer access.**

---

## Proposed Solution

### Phase 1: Define Stable Buffer Interface

Create a minimal, Eigen-free interface sufficient for Cython/NumPy:

```cpp
// Minimal interface required for NumPy compatibility
class MatrixBufferInterface {
public:
    double* data();
    const double* data() const;
    size_t rows() const;
    size_t cols() const;
    int innerStride() const;
    int outerStride() const;
    bool rowMajor() const;
    void resize(size_t rows, size_t cols);
    double getValue(size_t i, size_t j) const;
    void setValue(size_t i, size_t j, double value);
};
```

### Phase 2: Refactor Matrix<T> with PIMPL

**Public header (no Eigen):**

```cpp
// Matrix.h
#pragma once
#include <memory>
#include <cstddef>

namespace OpenMS {

namespace detail { template<typename T> struct MatrixImpl; }

template <typename Value>
class Matrix {
public:
    Matrix();
    Matrix(size_t rows, size_t cols, Value value = Value());
    ~Matrix();
    Matrix(const Matrix&);
    Matrix(Matrix&&) noexcept;
    Matrix& operator=(const Matrix&);
    Matrix& operator=(Matrix&&) noexcept;

    // Stable public API
    size_t rows() const;
    size_t cols() const;
    void resize(size_t rows, size_t cols);
    Value getValue(size_t i, size_t j) const;
    void setValue(size_t i, size_t j, Value value);

    // Buffer access for NumPy interop
    Value* data();
    const Value* data() const;
    int innerStride() const;
    int outerStride() const;
    bool rowMajor() const;

private:
    std::unique_ptr<detail::MatrixImpl<Value>> impl_;

    // Internal access for OpenMS code only
    friend class MatrixEigenAccess;
    void* eigenMatrix();
    const void* eigenMatrix() const;
};

extern template class Matrix<double>;
extern template class Matrix<float>;

} // namespace OpenMS
```

**Internal header (with Eigen, not installed):**

```cpp
// MatrixEigen.h - Internal use only
#pragma once
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <Eigen/Dense>

namespace OpenMS {

template <typename Value>
using EigenMatrixType = Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic>;

// Type-safe Eigen access for internal code
template <typename Value>
EigenMatrixType<Value>& getEigenMatrix(Matrix<Value>& mat) {
    return *static_cast<EigenMatrixType<Value>*>(mat.eigenMatrix());
}

template <typename Value>
const EigenMatrixType<Value>& getEigenMatrix(const Matrix<Value>& mat) {
    return *static_cast<const EigenMatrixType<Value>*>(mat.eigenMatrix());
}

// Zero-copy Eigen::Map view
template <typename Value>
auto eigenView(Matrix<Value>& mat) {
    return Eigen::Map<EigenMatrixType<Value>>(mat.data(), mat.rows(), mat.cols());
}

template <typename Value>
auto eigenView(const Matrix<Value>& mat) {
    return Eigen::Map<const EigenMatrixType<Value>>(mat.data(), mat.rows(), mat.cols());
}

} // namespace OpenMS
```

**Implementation (with Eigen):**

```cpp
// Matrix.cpp
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <Eigen/Dense>

namespace OpenMS::detail {

template <typename Value>
struct MatrixImpl {
    Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic> eigen_mat;

    MatrixImpl() = default;
    MatrixImpl(size_t rows, size_t cols, Value val) : eigen_mat(rows, cols) {
        eigen_mat.fill(val);
    }
};

} // namespace OpenMS::detail

namespace OpenMS {

template <typename Value>
Matrix<Value>::Matrix() : impl_(std::make_unique<detail::MatrixImpl<Value>>()) {}

template <typename Value>
Matrix<Value>::Matrix(size_t rows, size_t cols, Value value)
    : impl_(std::make_unique<detail::MatrixImpl<Value>>(rows, cols, value)) {}

template <typename Value>
Matrix<Value>::~Matrix() = default;

template <typename Value>
Value* Matrix<Value>::data() { return impl_->eigen_mat.data(); }

template <typename Value>
size_t Matrix<Value>::rows() const { return impl_->eigen_mat.rows(); }

template <typename Value>
size_t Matrix<Value>::cols() const { return impl_->eigen_mat.cols(); }

template <typename Value>
void* Matrix<Value>::eigenMatrix() { return &impl_->eigen_mat; }

// ... remaining implementations ...

template class Matrix<double>;
template class Matrix<float>;

} // namespace OpenMS
```

### Phase 3: Refactor TraceFitter and Similar Classes

Replace Eigen types in public signatures with standard types:

**Before:**
```cpp
class GenericFunctor {
    virtual int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) = 0;
    virtual int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) = 0;
};
```

**After:**
```cpp
class GenericFunctor {
    virtual int operator()(const std::vector<double>& x,
                          std::vector<double>& fvec) = 0;
    virtual int df(const std::vector<double>& x,
                  Matrix<double>& J) = 0;
};
```

Internal implementation adapts to Eigen:

```cpp
// TraceFitter.cpp
#include <unsupported/Eigen/NonLinearOptimization>

void TraceFitter::optimize_(std::vector<double>& x_init, GenericFunctor& functor) {
    // Wrap std::vector in Eigen::Map
    Eigen::Map<Eigen::VectorXd> x(x_init.data(), x_init.size());

    // Use internal adapter for LM solver
    // ... Eigen operations ...
}
```

### Phase 4: Handle BinnedSpectrum

Replace `Eigen::SparseVector` exposure with an opaque type:

```cpp
// BinnedSpectrum.h
class OPENMS_DLLAPI BinnedSpectrum {
public:
    // Instead of exposing SparseVectorType directly
    float getBinIntensity(size_t index) const;
    void setBinIntensity(size_t index, float value);
    std::vector<std::pair<size_t, float>> getNonZeroBins() const;
    size_t numNonZeroBins() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;  // Contains Eigen::SparseVector
};
```

---

## Migration Checklist

### Phase 1: Preparation
- [ ] Audit all public headers for Eigen includes/types
- [ ] Create `MatrixEigen.h` internal header
- [ ] Add CMake rules to not install internal headers

### Phase 2: Matrix Refactoring
- [ ] Implement PIMPL for `Matrix<T>`
- [ ] Add explicit template instantiations
- [ ] Update Cython `.pxd` declarations (minimal changes expected)
- [ ] Verify NumPy interop still works

### Phase 3: Algorithm Classes
- [ ] Refactor `TraceFitter::GenericFunctor`
- [ ] Refactor `LevMarqFitter1D`
- [ ] Refactor `NonNegativeLeastSquaresSolver`
- [ ] Update derived classes (`EGHTraceFitter`, `GaussTraceFitter`, `EmgFitter1D`)

### Phase 4: Sparse Types
- [ ] Refactor `BinnedSpectrum` to hide `SparseVector`
- [ ] Update any classes using `BinnedSpectrum::SparseVectorType`

### Phase 5: Cleanup
- [ ] Remove Eigen forward declarations from public headers
- [ ] Update documentation
- [ ] Add deprecation notices for any removed public APIs

---

## Cython Changes Required

The Cython bindings should require minimal changes since they already use the buffer interface:

```cython
# Matrix.pxd - No changes needed to declarations
cdef extern from "<OpenMS/DATASTRUCTURES/Matrix.h>" namespace "OpenMS":
    cdef cppclass Matrix[ValueT]:
        Matrix() except + nogil
        ValueT getValue(size_t i, size_t j) nogil
        void setValue(size_t i, size_t j, ValueT value) nogil
        size_t rows() nogil
        size_t cols() nogil
        int innerStride() nogil
        int outerStride() nogil
        bool rowMajor() nogil
        void resize(size_t rows, size_t cols) nogil
        ValueT* data() nogil
```

```cython
# MatrixDouble.pyx - No changes needed
def get_matrix_as_view(self):
    cdef _Matrix[double] * mat_ = self.inst.get()
    cdef double* data = mat_.data()
    # ... existing stride tricks code works unchanged ...
```

---

## Performance Considerations

| Aspect | Current | After Refactor |
|--------|---------|----------------|
| Matrix element access | Direct | +1 pointer indirection |
| Eigen operations (internal) | Optimal | Optimal (via `eigenView()`) |
| NumPy conversion | Zero-copy | Zero-copy (unchanged) |
| Binary size | Larger (Eigen in headers) | Smaller |
| Compile time | Slower (Eigen templates) | Faster for users |

**Key optimization:** Use `Eigen::Map` for zero-copy views:

```cpp
void compute(Matrix<double>& mat) {
    auto view = eigenView(mat);  // Zero-copy Eigen access
    view = view.array().exp();   // Full Eigen optimization
}
```

---

## Breaking Changes

### API Changes
1. `Matrix<T>` no longer inherits from `Eigen::Matrix`
2. `Matrix<T>::EigenMatrixType` typedef removed from public API
3. `Matrix<T>::getEigenMatrix()` returns opaque pointer (internal use only)
4. `TraceFitter::GenericFunctor` uses `std::vector` instead of `Eigen::VectorXd`
5. `BinnedSpectrum::SparseVectorType` no longer publicly accessible

### Migration for Downstream Users
```cpp
// Before
Matrix<double> mat(10, 10);
Eigen::MatrixXd& eigen = mat.getEigenMatrix();
eigen = eigen.transpose();

// After (if they need Eigen access)
#include <Eigen/Dense>
Matrix<double> mat(10, 10);
auto view = Eigen::Map<Eigen::MatrixXd>(mat.data(), mat.rows(), mat.cols());
view = view.transpose();
```

---

## Alternative Approaches Considered

### A. Keep Inheritance, Hide via Visibility
- Use `__attribute__((visibility("hidden")))` on Eigen symbols
- **Rejected:** Doesn't work reliably with templates

### B. Type Erasure
- Use `std::any` or `void*` for all matrix types
- **Rejected:** Too much runtime overhead, loses type safety

### C. Complete Eigen Removal
- Replace with custom matrix implementation
- **Rejected:** Loses Eigen's optimizations (SIMD, expression templates)

### D. Eigen::Map Only (No Storage)
- Store raw `std::vector`, use `Eigen::Map` everywhere
- **Viable alternative:** Slightly more code changes but no PIMPL overhead

---

## Timeline Estimate

| Phase | Effort | Dependencies |
|-------|--------|--------------|
| Phase 1: Preparation | 1-2 days | None |
| Phase 2: Matrix refactor | 3-5 days | Phase 1 |
| Phase 3: Algorithm classes | 3-5 days | Phase 2 |
| Phase 4: Sparse types | 2-3 days | Phase 2 |
| Phase 5: Cleanup & docs | 1-2 days | Phases 3-4 |

**Total: 10-17 days**

---

## References

- Current Matrix implementation: `src/openms/include/OpenMS/DATASTRUCTURES/Matrix.h`
- Cython bindings: `src/pyOpenMS/addons/MatrixDouble.pyx`
- TraceFitter: `src/openms/include/OpenMS/FEATUREFINDER/TraceFitter.h`
- BinnedSpectrum: `src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h`
- CMake visibility settings: `cmake/add_library_macros.cmake:165-166`
