#pragma once

#include <nanobind/nanobind.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
// IMSDataConsumer only forward-declares these, but the hand-off below needs them complete.
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <type_traits>

namespace nb = nanobind;

namespace
{
  /// Hands a Python callback an object that *owns* the value it is given.
  ///
  /// Spectra and chromatograms passed to a consumer live in the parser's batch buffer,
  /// which is cleared and refilled as soon as the callback returns. Handing Python a
  /// reference into that buffer therefore leaves any retained wrapper -- and every child
  /// alias derived from it, including the zero-copy ndarray from peaks_struct() --
  /// pointing at storage that has since been reused.
  ///
  /// Copying into the callback is not an option: IMSDataConsumer explicitly permits a
  /// consumer to modify the spectra, and MzMLHandler passes the possibly-modified value
  /// onward, so a copy would silently discard the consumer's edits.
  ///
  /// Moving instead removes the aliasing altogether. The value is moved out of the
  /// caller's storage for the duration of the call and moved back afterwards, so
  /// modification still works; if the consumer kept a reference, the value is copied back
  /// and Python is left holding a valid, independent snapshot instead of a dangling alias.
  /// MSSpectrum and MSChromatogram both have defaulted noexcept moves, so the normal path
  /// copies no peak data.
  template <typename T>
  class OwningCallbackArg
  {
    // The exception safety of restore() rests on these: the move-back cannot fail, so the
    // only throwing step is the copy construction on the retained path, which happens
    // before anything is assigned and therefore leaves the caller's storage untouched.
    static_assert(std::is_nothrow_move_assignable_v<T>,
                  "restore() relies on a non-throwing move assignment");
    static_assert(std::is_nothrow_move_constructible_v<T>,
                  "handing the value to Python relies on a non-throwing move construction");

  public:
    /// @param storage      the caller's value, moved out for the duration of the call.
    /// @param must_restore whether the value has to be put back afterwards. False only for
    ///        callers that provably discard it; see NanobindMSDataConsumer.
    OwningCallbackArg(T& storage, bool must_restore)
      : storage_(storage), obj_(nb::cast(std::move(storage), nb::rv_policy::move)),
        must_restore_(must_restore)
    {
      // nb::cast consults nanobind's instance registry for every policy except `copy`
      // (nb_type_put: `if (rvp != rv_policy::copy) { ... return the registered object; }`),
      // so a `move` cast returns an already-registered wrapper -- having moved nothing --
      // whenever the address is one Python already knows. No caller below can trigger that:
      // they all pass parser-owned batch storage that is never handed to Python. But if one
      // ever did, this class would quietly become a no-op and reinstate the very aliasing it
      // exists to remove, so say so instead. Nothing has moved in that case, which is why it
      // is safe to throw here: `storage` still holds its value.
      if (nb::inst_ptr<T>(obj_) == &storage)
      {
        throw std::runtime_error(
          "OwningCallbackArg: value is already registered with nanobind, so the hand-off "
          "would alias the caller's storage instead of owning a copy");
      }
    }

    OwningCallbackArg(const OwningCallbackArg&) = delete;
    OwningCallbackArg& operator=(const OwningCallbackArg&) = delete;

    /// Puts the possibly-modified value back into the caller's storage.
    ///
    /// Call this explicitly once the callback has returned normally. It is allowed to
    /// throw: if the value cannot be put back, the parser must not carry on with
    /// moved-from storage and pass that downstream, so the failure has to propagate.
    ///
    /// Does nothing when the caller never looks at the value again. That case is worth
    /// singling out because putting the value back costs a deep copy of the peak data
    /// whenever the consumer kept a reference -- straight into storage nobody reads.
    void restore()
    {
      if (restored_ || !must_restore_) { return; }
      T& edited = *nb::inst_ptr<T>(obj_); // a bound instance we just created; no conversion
      if (isUniquelyReferenced_())
      {
        storage_ = std::move(edited);
      }
      else
      {
        T copy(edited);            // copy first, so a throw here leaves storage_ untouched
        storage_ = std::move(copy);
      }
      restored_ = true;
    }

    /// Best-effort fallback for the exception path, where the callback's own error is
    /// already propagating and must not be masked by this one.
    ~OwningCallbackArg()
    {
      try
      {
        restore();
      }
      catch (...)
      {
      }
    }

    nb::object& handle() { return obj_; }

  private:
    bool isUniquelyReferenced_() const
    {
#if defined(NB_FREE_THREADED)
      // gil_scoped_acquire does not serialise threads on a free-threaded build, so the
      // reference count is not a sound uniqueness test there. Always take the copy path.
      return false;
#else
      return Py_REFCNT(obj_.ptr()) == 1;
#endif
    }

    T& storage_;
    nb::object obj_;
    bool must_restore_;
    bool restored_ = false;
  };
} // namespace

/// Duck-typing wrapper: bridges a Python consumer object (with consumeSpectrum,
/// consumeChromatogram, setExpectedSize, setExperimentalSettings methods) to the
/// C++ IMSDataConsumer interface. This enables MzMLFile.transform() / MzXMLFile.transform()
/// to call back into Python.
class NanobindMSDataConsumer : public OpenMS::Interfaces::IMSDataConsumer {
    nb::object py_consumer_;
    bool caller_reads_back_;
public:
    /// @param consumer the Python object to bridge.
    /// @param caller_reads_back whether the caller still looks at the spectrum or
    ///        chromatogram once the callback has returned. When it does, the value is put
    ///        back after every call, which costs a deep copy of the peak data whenever the
    ///        consumer kept a reference. The streaming overloads
    ///        MzMLFile/MzXMLFile::transform(filename, consumer) do not: they parse into a
    ///        function-local `PeakMap dummy` that is destroyed on the spot, and the handler
    ///        appends to it only when AlwaysAppendData is set (MzMLHandler.cpp:181-190,
    ///        MzXMLHandler.cpp:1251-1266), so nothing observable can read the value back and
    ///        the copy would be pure waste -- exactly for the `kept.append(s)` pattern that
    ///        motivates owning the value in the first place. Defaults to true, so a caller
    ///        added later is correct before it is fast.
    explicit NanobindMSDataConsumer(nb::object consumer, bool caller_reads_back = true)
      : py_consumer_(std::move(consumer)), caller_reads_back_(caller_reads_back) {}

    void consumeSpectrum(SpectrumType& s) override {
        nb::gil_scoped_acquire gil;
        try {
            // scoped inside the try: unwinding destroys it before the catch runs, while
            // nb::python_error still owns (and has cleared) the Python error indicator
            OwningCallbackArg<SpectrumType> arg(s, caller_reads_back_);
            py_consumer_.attr("consumeSpectrum")(arg.handle());
            arg.restore();
        } catch (nb::python_error& e) {
            e.restore();
            throw std::runtime_error("Python callback consumeSpectrum raised an exception");
        }
    }

    void consumeChromatogram(ChromatogramType& c) override {
        nb::gil_scoped_acquire gil;
        try {
            // scoped inside the try: unwinding destroys it before the catch runs, while
            // nb::python_error still owns (and has cleared) the Python error indicator
            OwningCallbackArg<ChromatogramType> arg(c, caller_reads_back_);
            py_consumer_.attr("consumeChromatogram")(arg.handle());
            arg.restore();
        } catch (nb::python_error& e) {
            e.restore();
            throw std::runtime_error("Python callback consumeChromatogram raised an exception");
        }
    }

    void setExpectedSize(size_t expectedSpectra, size_t expectedChromatograms) override {
        nb::gil_scoped_acquire gil;
        try {
            py_consumer_.attr("setExpectedSize")(expectedSpectra, expectedChromatograms);
        } catch (nb::python_error& e) {
            e.restore();
            throw std::runtime_error("Python callback setExpectedSize raised an exception");
        }
    }

    void setExperimentalSettings(const OpenMS::ExperimentalSettings& exp) override {
        nb::gil_scoped_acquire gil;
        try {
            py_consumer_.attr("setExperimentalSettings")(nb::cast(exp, nb::rv_policy::copy));
        } catch (nb::python_error& e) {
            e.restore();
            throw std::runtime_error("Python callback setExperimentalSettings raised an exception");
        }
    }
};
