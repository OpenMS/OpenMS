// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/Exception.h> // Include for exceptions

#include <vector>
#include <string>
#include <numeric>
#include <algorithm>

namespace OpenMS
{
  /**
    @brief Provides a non-owning view/adapter for an MSSpectrum containing centroided Ion Mobility data.

    This class wraps a pointer to an MSSpectrum (MSSpectrum*) that is assumed to be
    centroided and MUST contain a FloatDataArray with ion mobility data
    (as determined by MSSpectrum::containsIMData()) with drift time values
    corresponding to each peak.

    It provides convenient accessors for m/z, intensity, and drift time. Sorting
    operations modify the underlying MSSpectrum object directly.

    @warning This class stores a raw pointer (MSSpectrum*). The user is responsible
             for ensuring that the pointed-to MSSpectrum object remains valid for
             the lifetime of the IMFrameView object. Using a dangling pointer will
             lead to undefined behavior (likely crashes).

    @ingroup Kernel
  */
  class OPENMS_DLLAPI IMFrameView
  {
  public:
    /// Struct to hold peak data (m/z, intensity, drift time) for convenience
    /// (Identical to the previous version)
    struct IMPeak
    {
        double mz = 0.0;
        float intensity = 0.0f;
        double drift_time = 0.0; // Use double for external interface consistency

        IMPeak() = default;
        IMPeak(double mz, float intensity, double dt) :
          mz(mz), intensity(intensity), drift_time(dt) {}

        bool operator==(const IMPeak& other) const {
            return mz == other.mz && intensity == other.intensity && drift_time == other.drift_time;
        }
        bool operator!=(const IMPeak& other) const {
            return !(*this == other);
        }
    };

    ///@name Constructors and Destructor
    ///@{

    /// Constructor requires a pointer to an existing MSSpectrum.
    /// Performs validation checks on the pointed-to spectrum:
    /// - Pointer must not be null.
    /// - Spectrum must be centroided.
    /// - Spectrum must contain a FloatDataArray with ion mobility data (IMFormat::CONCATENATED).
    /// - The drift time array must have the same size as the number of peaks.
    /// Throws Exception::NullPointer, Exception::InvalidValue, or Exception::MissingInformation if validation fails.
    explicit IMFrameView(MSSpectrum* spectrum_ptr);

    // Default copy/move constructors/assignments are okay as they just copy/move the pointer.
    // The user remains responsible for the lifetime of the pointed-to object.
    IMFrameView(const IMFrameView& source) = default;
    IMFrameView(IMFrameView&& source) noexcept = default;
    IMFrameView& operator=(const IMFrameView& source) = default;
    IMFrameView& operator=(IMFrameView&& source) noexcept = default;


    /// Destructor (does not delete the pointed-to spectrum)
    virtual ~IMFrameView();
    ///@}

    ///@name Accessors (Const methods provide read-only access)
    ///@{
    /// Returns a const pointer to the underlying MSSpectrum.
    const MSSpectrum* getSpectrum() const;

    /// Returns a pointer to the underlying MSSpectrum.
    /// @warning Use caution when modifying the spectrum directly, as it might
    /// invalidate the assumptions of the IMFrameView (e.g., drift time consistency)
    /// if not done carefully. Prefer IMFrameView methods where possible.
    MSSpectrum* getSpectrum();

    /// Read-only access to the peaks as IMPeak structs. Constructs these on the fly. (const)
    std::vector<IMPeak> getPeaks() const;

    /// Read-only access to m/z values (const)
    std::vector<double> getMZArray() const;

    /// Read-only access to intensity values (const)
    std::vector<float> getIntensityArray() const;

    /// Read-only access to drift time values (converted to double) (const)
    std::vector<double> getDriftTimeArray() const;

    /// Returns the number of peaks in the frame (const)
    Size size() const;

    /// Returns true if the frame is empty (no peaks) (const)
    bool empty() const;

    /// Returns true if the object is valid (const)
    /// Note: The constructor already validates the spectrum, so this will always return true
    /// for a successfully constructed object.
    bool isValid() const; 

    /// Gets a const pointer to the FloatDataArray containing drift times. (const)
    /// Returns nullptr if not found or if name mismatch. Size check is separate.
    const MSSpectrum::FloatDataArray* getDriftTimeDataArrayPtr() const;

    /// Returns the retention time (const)
    double getRT() const;

    /// Sets the retention time on the underlying spectrum (non-const)
    void setRT(double rt);

    /// Returns the MS level (const)
    UInt getMSLevel() const;

    /// Sets the MS level on the underlying spectrum (non-const)
    void setMSLevel(UInt level);

    /// Gets the native ID string (const)
    const std::string& getNativeID() const;

    /// Sets the native ID string on the underlying spectrum (non-const)
    void setNativeID(const std::string& id);

    // Add other relevant getters (const) / setters (non-const) as needed
    ///@}

    ///@name Operations (Modify the underlying spectrum)
    ///@{

    /// Sorts peaks by m/z value (ascending). Also sorts drift times accordingly. (non-const)
    /// Modifies the underlying MSSpectrum and its FloatDataArray.
    void sortByMZ();

    /// Sorts peaks by intensity. Also sorts drift times accordingly. (non-const)
    /// Modifies the underlying MSSpectrum and its FloatDataArray.
    void sortByIntensity(bool reverse = true);

    /// Sorts peaks by drift time (ascending). Also sorts drift times accordingly. (non-const)
    /// Modifies the underlying MSSpectrum and its FloatDataArray.
    void sortByDriftTime();

    ///@}

    ///@name Comparison operators (Compare pointers for identity, or underlying spectra for equality)
    ///@{
    /// Equality operator: Compares the underlying spectra for value equality.
    bool operator==(const IMFrameView& rhs) const;

    /// Inequality operator
    bool operator!=(const IMFrameView& rhs) const;
    ///@}

  private:
    /// Pointer to the non-owned MSSpectrum object.
    MSSpectrum* spectrum_ptr_;

    /// Helper to find the drift time array (const version)
    const MSSpectrum::FloatDataArray* findDriftTimeArray_() const;

    /// Helper to find the drift time array (non-const version for sorting)
    MSSpectrum::FloatDataArray* findDriftTimeArray_();

    /// Performs validation checks on the pointed-to spectrum. Called by constructor. (const)
    void validateSpectrum_(const std::string& context_message = "") const;

    /// Internal sorting helper (operates on *spectrum_ptr_ and its drift time array) (non-const)
    template <typename Compare> 
    void sortPeaksAndDriftTime_(Compare comp)
    {
        checkPointer_();
        
        MSSpectrum* spec = spectrum_ptr_;
        
        // Create a vector of indices
        std::vector<Size> indices(spec->size());
        std::iota(indices.begin(), indices.end(), 0);
        
        // Sort the indices based on the comparison function
        std::stable_sort(indices.begin(), indices.end(), 
                        [&](Size i, Size j) { return comp(i, j); });
        
        // Use MSSpectrum::select to reorder the peaks and all data arrays
        spec->select(indices);
    }

    /// Internal helper to check pointer validity
    inline void checkPointer_() const
    {
        // This check should ideally only fail if the object was moved from
        // or if constructed with nullptr (which constructor prevents).
        // It doesn't protect against the underlying object being deleted elsewhere.
        if (!spectrum_ptr_) {
             throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Internal spectrum pointer is null. View is invalid (e.g., moved from?).");
        }
    }


  }; // class IMFrameView

} // namespace OpenMS
