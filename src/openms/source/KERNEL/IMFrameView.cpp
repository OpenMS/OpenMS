// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/IMFrameView.h>

namespace OpenMS
{
  IMFrameView::IMFrameView(MSSpectrum* spectrum_ptr) : spectrum_ptr_(spectrum_ptr)
  {
    validateSpectrum_("Constructor");
  }

  IMFrameView::~IMFrameView() = default;

  const MSSpectrum* IMFrameView::getSpectrum() const
  {
    checkPointer_();
    return spectrum_ptr_;
  }

  MSSpectrum* IMFrameView::getSpectrum()
  {
    checkPointer_();
    return spectrum_ptr_;
  }

  std::vector<IMFrameView::IMPeak> IMFrameView::getPeaks() const
  {
    checkPointer_();
    
    const MSSpectrum* spec = spectrum_ptr_;
    const MSSpectrum::FloatDataArray* drift_times = findDriftTimeArray_();
    
    std::vector<IMPeak> result;
    result.reserve(spec->size());
    
    for (Size i = 0; i < spec->size(); ++i)
    {
      const Peak1D& peak = (*spec)[i];
      result.emplace_back(peak.getMZ(), peak.getIntensity(), static_cast<double>((*drift_times)[i]));
    }
    
    return result;
  }

  std::vector<double> IMFrameView::getMZArray() const
  {
    checkPointer_();
    
    const MSSpectrum* spec = spectrum_ptr_;
    std::vector<double> result;
    result.reserve(spec->size());
    
    for (Size i = 0; i < spec->size(); ++i)
    {
      result.push_back((*spec)[i].getMZ());
    }
    
    return result;
  }

  std::vector<float> IMFrameView::getIntensityArray() const
  {
    checkPointer_();
    
    const MSSpectrum* spec = spectrum_ptr_;
    std::vector<float> result;
    result.reserve(spec->size());
    
    for (Size i = 0; i < spec->size(); ++i)
    {
      result.push_back((*spec)[i].getIntensity());
    }
    
    return result;
  }

  std::vector<double> IMFrameView::getDriftTimeArray() const
  {
    checkPointer_();
    
    const MSSpectrum::FloatDataArray* drift_times = findDriftTimeArray_();
    std::vector<double> result;
    result.reserve(drift_times->size());
    
    for (Size i = 0; i < drift_times->size(); ++i)
    {
      result.push_back(static_cast<double>((*drift_times)[i]));
    }
    
    return result;
  }

  Size IMFrameView::size() const
  {
    checkPointer_();
    return spectrum_ptr_->size();
  }

  bool IMFrameView::empty() const
  {
    checkPointer_();
    return spectrum_ptr_->empty();
  }

  bool IMFrameView::isValid() const
  {
    try
    {
      validateSpectrum_("isValid");
      return true;
    }
    catch (...)
    {
      return false;
    }
  }

  const MSSpectrum::FloatDataArray* IMFrameView::getDriftTimeDataArrayPtr() const
  {
    checkPointer_();
    return findDriftTimeArray_();
  }

  double IMFrameView::getRT() const
  {
    checkPointer_();
    return spectrum_ptr_->getRT();
  }

  void IMFrameView::setRT(double rt)
  {
    checkPointer_();
    spectrum_ptr_->setRT(rt);
  }

  UInt IMFrameView::getMSLevel() const
  {
    checkPointer_();
    return spectrum_ptr_->getMSLevel();
  }

  void IMFrameView::setMSLevel(UInt level)
  {
    checkPointer_();
    spectrum_ptr_->setMSLevel(level);
  }

  const std::string& IMFrameView::getNativeID() const
  {
    checkPointer_();
    return spectrum_ptr_->getNativeID();
  }

  void IMFrameView::setNativeID(const std::string& id)
  {
    checkPointer_();
    spectrum_ptr_->setNativeID(id);
  }

  void IMFrameView::sortByMZ()
  {
    checkPointer_();
    
    // Use the template method with a comparator for m/z
    sortPeaksAndDriftTime_([spec = spectrum_ptr_](Size i, Size j) {
      return (*spec)[i].getMZ() < (*spec)[j].getMZ();
    });
  }

  void IMFrameView::sortByIntensity(bool reverse)
  {
    checkPointer_();
    
    // Use the template method with a comparator for intensity
    if (reverse)
    {
      sortPeaksAndDriftTime_([spec = spectrum_ptr_](Size i, Size j) {
        return (*spec)[i].getIntensity() > (*spec)[j].getIntensity();
      });
    }
    else
    {
      sortPeaksAndDriftTime_([spec = spectrum_ptr_](Size i, Size j) {
        return (*spec)[i].getIntensity() < (*spec)[j].getIntensity();
      });
    }
  }

  void IMFrameView::sortByDriftTime()
  {
    checkPointer_();
    
    // Use the template method with a comparator for drift time
    MSSpectrum::FloatDataArray* dt_array = findDriftTimeArray_();
    sortPeaksAndDriftTime_([dt_array](Size i, Size j) {
      return (*dt_array)[i] < (*dt_array)[j];
    });
  }

  bool IMFrameView::operator==(const IMFrameView& rhs) const
  {
    // First check if they point to the same spectrum
    if (spectrum_ptr_ == rhs.spectrum_ptr_)
    {
      return true;
    }
    
    // If either pointer is null, they can't be equal (since they're not the same pointer)
    if (!spectrum_ptr_ || !rhs.spectrum_ptr_)
    {
      return false;
    }
    
    // Check if the spectra have the same number of peaks
    if (spectrum_ptr_->size() != rhs.spectrum_ptr_->size())
    {
      return false;
    }
    
    // Get the drift time arrays
    const MSSpectrum::FloatDataArray* dt_array1 = findDriftTimeArray_();
    const MSSpectrum::FloatDataArray* dt_array2 = rhs.findDriftTimeArray_();
    
    // Check if both have valid drift time arrays of the same size
    if (!dt_array1 || !dt_array2 || dt_array1->size() != dt_array2->size())
    {
      return false;
    }
    
    // Compare each peak and its drift time
    for (Size i = 0; i < spectrum_ptr_->size(); ++i)
    {
      const Peak1D& peak1 = (*spectrum_ptr_)[i];
      const Peak1D& peak2 = (*rhs.spectrum_ptr_)[i];
      
      // Compare m/z, intensity, and drift time
      if (peak1.getMZ() != peak2.getMZ() ||
          peak1.getIntensity() != peak2.getIntensity() ||
          (*dt_array1)[i] != (*dt_array2)[i])
      {
        return false;
      }
    }
    
    // All peaks and drift times are equal
    return true;
  }

  bool IMFrameView::operator!=(const IMFrameView& rhs) const
  {
    return !(*this == rhs);
  }

  const MSSpectrum::FloatDataArray* IMFrameView::findDriftTimeArray_() const
  {
    checkPointer_();
    
    const MSSpectrum* spec = spectrum_ptr_;
    try
    {
      // Get the ion mobility data array index
      Size im_index = spec->getIMData().first;
      return &(spec->getFloatDataArrays()[im_index]);
    }
    catch (Exception::MissingInformation&)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Could not find ion mobility data array");
    }
  }

  MSSpectrum::FloatDataArray* IMFrameView::findDriftTimeArray_()
  {
    checkPointer_();
    
    MSSpectrum* spec = spectrum_ptr_;
    try
    {
      // Get the ion mobility data array index
      Size im_index = spec->getIMData().first;
      return &(spec->getFloatDataArrays()[im_index]);
    }
    catch (Exception::MissingInformation&)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Could not find ion mobility data array");
    }
  }

  void IMFrameView::validateSpectrum_(const std::string& context_message) const
  {
    // Check if the pointer is null
    if (!spectrum_ptr_)
    {
      throw Exception::NullPointer(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    
    // Check if the spectrum is centroided
    if (spectrum_ptr_->getType() != SpectrumSettings::CENTROID)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                  "Spectrum must be centroided. " + context_message,
                                  "non-centroided");
    }
    
    // Check if the spectrum has the correct ion mobility format (CONCATENATED)
    if (IMTypes::determineIMFormat(*spectrum_ptr_) != IMFormat::CONCATENATED)
    {
      IMFormat format = IMTypes::determineIMFormat(*spectrum_ptr_);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Spectrum must have ion mobility data in CONCATENATED format (with drift time array). " + context_message,
                                  toString(format));
    }
    
    // Find the drift time array
    try
    {
      // Get the ion mobility data array
      Size im_index = spectrum_ptr_->getIMData().first;
      const MSSpectrum::FloatDataArray* dt_array = &(spectrum_ptr_->getFloatDataArrays()[im_index]);
      
      // Check if the drift time array has the same size as the number of peaks
      if (dt_array->size() != spectrum_ptr_->size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "Drift time array size (" + std::to_string(dt_array->size()) + 
                                     ") must match the number of peaks (" + std::to_string(spectrum_ptr_->size()) + 
                                     "). " + context_message,
                                     std::to_string(dt_array->size()));
      }
    }
    catch (Exception::MissingInformation&)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Spectrum must contain an ion mobility data array. " + context_message);
    }
  }


} // namespace OpenMS