// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Administrator $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/RangeManager.h>
#include <map>
#include <set>

namespace OpenMS
{
  class MSSpectrum;  // Forward declaration for MSSpectrum
  /**
    @brief Range manager for MS spectra with MS level-specific range management
    
    This class manages m/z, intensity, and ion mobility ranges for MS spectra,
    with the ability to filter ranges by MS level.
    
    @ingroup Kernel
  */
  class OPENMS_DLLAPI SpectrumRangeManager : public RangeManager<RangeMZ, RangeIntensity, RangeMobility, RangeRT>
  {
  public:
    /// Base type
    using BaseType = RangeManager<RangeMZ, RangeIntensity, RangeMobility, RangeRT>;
    
    /// Default constructor
    SpectrumRangeManager() = default;
    
    /// Copy constructor
    SpectrumRangeManager(const SpectrumRangeManager& source) = default;
    
    /// Move constructor
    SpectrumRangeManager(SpectrumRangeManager&& source) = default;
    
    /// Assignment operator
    SpectrumRangeManager& operator=(const SpectrumRangeManager& source) = default;
    
    /// Move assignment operator
    SpectrumRangeManager& operator=(SpectrumRangeManager&& source) = default;
    
    /// Destructor
    ~SpectrumRangeManager() = default;
    
    /**
      @brief Clears all ranges (global and MS level-specific)
    */
    void clear()
    {
      BaseType::clearRanges();
      ms_level_ranges_.clear();
    }
    
    /**
      @brief Extends the ranges with the ranges of another range manager
      
      @param other The other range manager to extend from
      @param ms_level The MS level for which to extend the ranges (0 for global ranges)
    */
    void extend(const BaseType& other, UInt ms_level = 0)
    {
      if (ms_level == 0)
      {
        BaseType::extendUnsafe(other);
      }
      else
      {
        if (ms_level_ranges_.find(ms_level) == ms_level_ranges_.end())
        {
          ms_level_ranges_[ms_level] = other;
        }
        else
        {
          ms_level_ranges_[ms_level].extendUnsafe(other);
        }
      }
    }
    
    /**
      @brief Gets the ranges for a specific MS level
      
      @param ms_level The MS level for which to retrieve the ranges (0 for global ranges)
      @return The ranges for the specified MS level (returns the global ranges if no specific ranges exist for the requested MS level)
    */
    const BaseType& operator()(UInt ms_level = 0) const
    {
      if (ms_level == 0 || ms_level_ranges_.find(ms_level) == ms_level_ranges_.end())
      {
        return *this;
      }
      return ms_level_ranges_.at(ms_level);
    }

    /**
      @brief Gets all MS levels for which specific ranges exist
      
      @return The set of MS levels
    */
    std::set<UInt> getMSLevels() const
    {
      std::set<UInt> ms_levels;
      for (const auto& [level, _] : ms_level_ranges_)
      {
        ms_levels.insert(level);
      }
      return ms_levels;
    }

    /**
      @brief Extends the RT range with an MS level parameter
      
      @param rt The RT value to extend with
      @param ms_level The MS level for which to extend the RT range (0 for global range)
    */
    void extendRT(double rt, UInt ms_level)
    {
      if (ms_level == 0)
      {
        BaseType::extendRT(rt);
      }
      else
      {
        if (ms_level_ranges_.find(ms_level) == ms_level_ranges_.end())
        {
          ms_level_ranges_[ms_level] = BaseType();
        }
        ms_level_ranges_[ms_level].extendRT(rt);
      }
    }

    /**
      @brief Extends the RT range with an MS level parameter
      
      @param rt The RT value to extend with
      @param ms_level The MS level for which to extend the RT range (0 for global range)
    */
    void extendMZ(double rt, UInt ms_level)
    {
      if (ms_level == 0)
      {
        BaseType::extendMZ(rt);
      }
      else
      {
        if (ms_level_ranges_.find(ms_level) == ms_level_ranges_.end())
        {
          ms_level_ranges_[ms_level] = BaseType();
        }
        ms_level_ranges_[ms_level].extendMZ(rt);
      }
    }


    /**
      @brief Extends the ranges with the ranges of a spectrum using an MS level parameter
      
      @param spectrum The spectrum whose ranges to extend from
      @param ms_level The MS level for which to extend the ranges (0 for global ranges)
    */
    void extendUnsafe(const MSSpectrum& spectrum, UInt ms_level)
    {
      if (ms_level == 0)
      {
        BaseType::extendUnsafe(spectrum);
      }
      else
      {
        if (ms_level_ranges_.find(ms_level) == ms_level_ranges_.end())
        {
          ms_level_ranges_[ms_level] = BaseType();
        }
        ms_level_ranges_[ms_level].extendUnsafe(spectrum);
      }
    }

  protected:
    /// MS level-specific ranges
    std::map<UInt, BaseType> ms_level_ranges_;
  };

} // namespace OpenMS