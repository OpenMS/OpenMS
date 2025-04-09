// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <functional>
#include <algorithm>
#include <vector>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  /**
    @defgroup RangeUtils RangeUtils

    @brief Predicates for range operations

    @ingroup Kernel
    A group of predicates that can be used to perform range operations on MS data.
    They operate on classes that have the save interface as Spectrum or Peak1D or Peak2D, respectively.
    <BR>
    <BR>
    The code for the removal of spectra in a certain retention time range from a vector of spectra might look like this:

    @code
    //data
    std::vector< MSSpectrum> spectra;

    //... spectra are added to the vector ...

    //range from 0.0 to 36.0 s
    InRTRange< MSSpectrum> range(0.0, 36.0);

    //remove the range
    spectra.erase(remove_if(spectra.begin(), spectra.end(), range), spectra.end());
    @endcode

    The code for the removal of peaks within certain intensity range from a spectrum might look like this:

    @code
    //data
    MSSpectrum spectrum;

    //... peaks are added to the spectrum ...

    //range from 0.0 to 5000.0 intensity
    InIntensityRange range< Peak1D >(0.0, 5000.0);

    //remove the range
    spectrum.erase(remove_if(spectrum.begin(), spectrum.end(), range), spectrum.end());
    @endcode
  */

  /**
    @brief Predicate that determines if a class has a certain metavalue

    MetaContainer must be a MetaInfoInterface or have the same interface

    @ingroup RangeUtils
  */
  template <class MetaContainer>
  class HasMetaValue
  {
public:
    /**
      @brief Constructor

      @param metavalue MetaValue that needs to be present.
      @param reverse if @p reverse is true, operator() returns true if the metavalue does not exist.
    */
    HasMetaValue(String metavalue, bool reverse = false) :
      metavalue_key_(metavalue),
      reverse_(reverse)
    {}

    inline bool operator()(const MetaContainer& s) const
    {
      bool has_meta_value = s.metaValueExists(metavalue_key_);
      // XOR(^): same as 'if (rev_) return !(test) else return test;' where (test) is the condition;   Speed: XOR is about 25% faster in VS10
      return reverse_ ^ has_meta_value;
    }

protected:
    String metavalue_key_;
    bool reverse_;
  };


  /**
    @brief Predicate that determines if a spectrum lies inside/outside a specific retention time range

    SpectrumType must be a Spectrum or have the same interface

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class InRTRange
  {
public:
    /**
      @brief Constructor

      @param min lower boundary
      @param max upper boundary
      @param reverse if @p reverse is true, operator() returns true if the spectrum lies outside the
      range
    */
    InRTRange(double min, double max, bool reverse = false) :
      min_(min),
      max_(max),
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      double tmp = s.getRT();
      // XOR(^): same as 'if (rev_) return !(test) else return test;' where (test) is the condition;   Speed: XOR is about 25% faster in VS10
      return reverse_ ^ (min_ <= tmp && tmp <= max_);
    }

protected:
    double min_, max_;
    bool reverse_;
  };

  /**
    @brief Predicate that determines if a spectrum lies inside/outside a specific MS level set

    SpectrumType must be a Spectrum or have the same interface

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class InMSLevelRange
  {
public:
    /**
      @brief Constructor

      @param levels an array of MS levels
      @param reverse if @p reverse is true, operator() returns true if the spectrum lies outside the
      set
    */
    InMSLevelRange(const IntList& levels, bool reverse = false) :
      levels_(levels),
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      Int tmp = s.getMSLevel();
      // XOR(^): same as 'if (rev_) return !(test) else return test;' where (test) is the condition;   Speed: XOR is about 25% faster in VS10
      return reverse_ ^ (std::find(levels_.begin(), levels_.end(), tmp) != levels_.end());
    }

protected:
    IntList levels_;
    bool reverse_;
  };

  /**
    @brief Predicate that determines if a spectrum has a certain scan mode

    SpectrumType must be a Spectrum or have the same interface (SpectrumSettings)

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class HasScanMode
  {
public:
    /**
      @brief Constructor

      @param mode scan mode
      @param reverse if @p reverse is true, operator() returns true if the spectrum has a different
      scan mode
    */
    HasScanMode(Int mode, bool reverse = false) :
      mode_(mode),
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      // XOR(^): same as 'if (rev_) return !(test) else return test;' where (test) is the condition;   Speed: XOR is about 25% faster in VS10
      return reverse_ ^ (s.getInstrumentSettings().getScanMode() == mode_);
    }

protected:
    Int mode_;
    bool reverse_;
  };

  /**
    @brief Predicate that determines if a spectrum has a certain scan polarity

    SpectrumType must be a Spectrum or have the same interface (SpectrumSettings)

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class HasScanPolarity
  {
public:
    /**
      @brief Constructor

      @param polarity scan polarity
      @param reverse if @p reverse is true, operator() returns true if the spectrum has a different
      scan polarity
    */
    HasScanPolarity(Int polarity, bool reverse = false) :
      polarity_(polarity),
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      // XOR(^): same as 'if (rev_) return !(test) else return test;' where (test) is the condition;   Speed: XOR is about 25% faster in VS10
      return reverse_ ^ (s.getInstrumentSettings().getPolarity() == polarity_);
    }

protected:
    Int polarity_;
    bool reverse_;
  };


  /**
    @brief Predicate that determines if a spectrum is empty.

    SpectrumType must have a empty() method

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class IsEmptySpectrum
  {
public:
    /**
      @brief Constructor

      @param reverse if @p reverse is true, operator() returns true if the spectrum is not empty
    */
    explicit IsEmptySpectrum(bool reverse = false) :
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      // XOR(^): same as 'if (rev_) return !(test) else return test;' where (test) is the condition;   Speed: XOR is about 25% faster in VS10
      return reverse_ ^ s.empty();
    }

protected:
    bool reverse_;
  };

  /**
    @brief Predicate that determines if a spectrum is a zoom (enhanced resolution) spectrum.

    SpectrumType must have a getInstrumentSettings() method

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class IsZoomSpectrum
  {
public:
    /**
      @brief Constructor

      @param reverse if @p reverse is true, operator() returns true if the spectrum is not a zoom
      spectrum
    */
    explicit IsZoomSpectrum(bool reverse = false) :
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      // XOR(^): same as 'if (rev_) return !(test) else return test;' where (test) is the condition;   Speed: XOR is about 25% faster in VS10
      return reverse_ ^ s.getInstrumentSettings().getZoomScan();
    }

protected:
    bool reverse_;
  };


  /**
    @brief Predicate that determines if a spectrum was generated using any activation method given
    in the constructor list.

    SpectrumType must have a getPrecursors() method

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class HasActivationMethod
  {
public:
    /**
      @brief Constructor

      @param methods List of methods that is compared against precursor activation methods.
      @param reverse if @p reverse is true, operator() returns true if the spectrum is not using one
      of the specified activation methods.
    */
    HasActivationMethod(const StringList& methods, bool reverse = false) :
      methods_(methods),
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      for (const Precursor& p : s.getPrecursors())
      {
        for (const Precursor::ActivationMethod am : p.getActivationMethods())
        {
          if (ListUtils::contains(methods_, Precursor::NamesOfActivationMethod[am]))
          {
            // found matching activation method
            if (reverse_) return false;
            else return true;
          }
        }
      }

      return reverse_;
    }

protected:
    StringList methods_;
    bool reverse_;
  };

  /**
    @brief Predicate that determines if a spectrum's precursor is within a certain m/z range.

    SpectrumType must have a getPrecursors() method

    If multiple precursors are present, all must fulfill the range criterion.

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class InPrecursorMZRange
  {
public:
    /**
      @brief Constructor

      @param mz_left left m/z boundary (closed interval)
      @param mz_right right m/z boundary (closed interval)
      @param reverse if @p reverse is true, operator() returns true if the precursor's m/z is outside of the given interval.
    */
    InPrecursorMZRange(const double& mz_left, const double& mz_right, bool reverse = false) :
      mz_left_(mz_left),
      mz_right_(mz_right),
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      for (std::vector<Precursor>::const_iterator it = s.getPrecursors().begin(); it != s.getPrecursors().end(); ++it)
      {
        //std::cerr << mz_left_ << " " << mz_right_ << " " << it->getMZ() << "\n";
        if (!(mz_left_ <= it->getMZ() && it->getMZ() <= mz_right_))
        { // found PC outside of allowed window
          if (reverse_) return true;
          else return false;
        }
      }

      if (reverse_) return false;
      else return true;
    }

protected:
    double mz_left_;
    double mz_right_;
    bool reverse_;
  };


  /**
    @brief Predicate that determines if a spectrum has a certain precursor charge as given in the
    constructor list.

    SpectrumType must have a getPrecursors() method

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class HasPrecursorCharge
  {
public:
    /**
      @brief Constructor

      @param charges List of charges that is compared against precursor charge.
      @param reverse if @p reverse is true, operator() returns true if the spectrum has not one of
      the specified precursor charges.
    */
    HasPrecursorCharge(const IntList& charges, bool reverse = false) :
      charges_(charges),
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      bool match = false;
      for (std::vector<Precursor>::const_iterator it = s.getPrecursors().begin(); it != s.getPrecursors().end(); ++it)
      {
        Int tmp = it->getCharge();
        match = match || (std::find(charges_.begin(), charges_.end(), tmp) != charges_.end());
      }

      if (reverse_) return !match;
      else return match;
    }

protected:
    IntList charges_;
    bool reverse_;
  };


  /**
    @brief Predicate that determines if a peak lies inside/outside a specific m/z range

    PeakType must have a getPosition() method.

    @note It is assumed that the m/z dimension is dimension 0!

    @ingroup RangeUtils
  */
  template <class PeakType>
  class InMzRange
  {
public:
    /**
      @brief Constructor

      @param min lower boundary
      @param max upper boundary
      @param reverse if @p reverse is true, operator() returns true if the peak lies outside the
      range
    */
    InMzRange(double min, double max, bool reverse = false) :
      min_(min),
      max_(max),
      reverse_(reverse)
    {}

    inline bool operator()(const PeakType& p) const
    {
      double tmp = p.getPosition()[0];
      // XOR(^): same as 'if (rev_) return !(test) else return test;' where (test) is the condition;   Speed: XOR is about 25% faster in VS10
      return reverse_ ^ (min_ <= tmp && tmp <= max_);
    }

protected:
    double min_, max_;
    bool reverse_;
  };

  /**
    @brief Predicate that determines if a peak lies inside/outside a specific intensity range

    PeakType must have a getIntensity() method.

    @ingroup RangeUtils
  */
  template <class PeakType>
  class InIntensityRange
  {
public:
    /**
      @brief Constructor

      @param min lower boundary
      @param max upper boundary
      @param reverse if @p reverse is true, operator() returns true if the peak lies outside the set
    */
    InIntensityRange(double min, double max, bool reverse = false) :
      min_(min),
      max_(max),
      reverse_(reverse)
    {}

    inline bool operator()(const PeakType& p) const
    {
      double tmp = p.getIntensity();
      // XOR(^): same as 'if (rev_) return !(test) else return test;' where (test) is the condition;   Speed: XOR is about 25% faster in VS10
      return reverse_ ^ (min_ <= tmp && tmp <= max_);
    }

protected:
    double min_, max_;
    bool reverse_;
  };

  /**
    @brief Predicate that determines if an MSn spectrum was generated with a collision energy in the given range.
    @note This applies only to CID and HCD spectra. For spectra that do not have a collision energy, the predicate will return true.
    @note This predicate will return always true for spectra with getMSLevel() = 1.

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class IsInCollisionEnergyRange
  {
public:
    /**
      @brief Constructor

      @param min minimum collision energy to be included in the range.
      @param max maximum collision energy to be included in the range.
      @param reverse if @p reverse is true, operator() returns true if the collision energy lies outside the range.
    */
    IsInCollisionEnergyRange(double min, double max, bool reverse = false) :
      min_energy_(min),
      max_energy_(max),
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      // leave non-fragmentation spectra untouched
      if (s.getMSLevel() == 1) return false;

      bool isIn = false;
      bool hasCollisionEnergy = false;
      for (std::vector<Precursor>::const_iterator it = s.getPrecursors().begin(); it != s.getPrecursors().end(); ++it)
      {
        if (it->metaValueExists("collision energy"))
        {
          hasCollisionEnergy = true;
          double cE = it->getMetaValue("collision energy");
          isIn |= !(cE > max_energy_ || cE < min_energy_);
        }
      }

      // we accept all spectra that have no collision energy value
      if (!hasCollisionEnergy) return false;

      if (reverse_) return !isIn;
      else return isIn;
    }

private:
    double min_energy_, max_energy_;
    bool reverse_;
  };

  /**
    @brief Predicate that determines if the width of the isolation window of an MSn spectrum is in the given range.
    @note This predicate will return always true for spectra with getMSLevel() = 1.

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class IsInIsolationWindowSizeRange
  {

public:
    /**
      @brief Constructor

      @param min_size minimum width of the isolation window.
      @param max_size maximum width of the isolation window.
      @param reverse if @p reverse is true, operator() returns true if the width of the isolation window lies outside the range.
    */
    IsInIsolationWindowSizeRange(double min_size, double max_size, bool reverse = false) :
      min_size_(min_size),
      max_size_(max_size),
      reverse_(reverse)
    {}

    inline bool operator()(const SpectrumType& s) const
    {
      // leave non-fragmentation spectra untouched
      if (s.getMSLevel() == 1) return false;

      bool isIn = false;
      for (std::vector<Precursor>::const_iterator it = s.getPrecursors().begin(); it != s.getPrecursors().end(); ++it)
      {
        const double isolationWindowSize = it->getIsolationWindowUpperOffset() + it->getIsolationWindowLowerOffset();
        isIn |= !(isolationWindowSize > max_size_ || isolationWindowSize < min_size_);
      }

      if (reverse_) return !isIn;
      else return isIn;
    }

private:
    double min_size_, max_size_;
    bool reverse_;
  };

    /**
    @brief Predicate that determines if the isolation window covers ANY of the given m/z values.
    @note This predicate will return always true for spectra with getMSLevel() = 1.

    @ingroup RangeUtils
  */
  template <class SpectrumType>
  class IsInIsolationWindow
  {

public:
    /**
      @brief Constructor

      @param vec_mz Vector of m/z values, of which at least one needs to be covered
      @param reverse if @p reverse is true, operator() returns true if the isolation window is outside for ALL m/z values.
    */
    IsInIsolationWindow(std::vector<double> vec_mz, bool reverse = false) :
      vec_mz_(vec_mz),
      reverse_(reverse)
    {
      std::sort(vec_mz_.begin(), vec_mz_.end());
    }

    inline bool operator()(const SpectrumType& s) const
    {
      // leave non-fragmentation spectra untouched
      if (s.getMSLevel() == 1) return false;

      bool isIn = false;
      for (std::vector<Precursor>::const_iterator it = s.getPrecursors().begin(); it != s.getPrecursors().end(); ++it)
      {
        if (it->getIsolationWindowLowerOffset() == 0 || it->getIsolationWindowUpperOffset() == 0)
        {
          OPENMS_LOG_WARN << "IsInIsolationWindow(): Lower/Upper Offset for Precursor Isolation Window is Zero! " << 
            "Filtering will probably be too strict (unless you hit the exact precursor m/z)!" << std::endl;
        }
        const double lower_mz = it->getMZ() - it->getIsolationWindowLowerOffset();
        std::vector<double>::const_iterator it_mz = std::lower_bound(vec_mz_.begin(), vec_mz_.end(), lower_mz);
        if (it_mz != vec_mz_.end()) // left side ok
        { // right side?
          const double upper_mz = it->getMZ() + it->getIsolationWindowUpperOffset();
          isIn |= (*it_mz <= upper_mz);
        }
      }

      if (reverse_) return !isIn;
      else return isIn;
    }

private:
    std::vector<double> vec_mz_;
    bool reverse_;
  };

} // namespace OpenMS

