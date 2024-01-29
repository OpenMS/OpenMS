// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/DPeak.h>

namespace OpenMS
{

  /**
    @brief Abstract base class for all D-dimensional models.

    Every derived class has to implement the static functions
    "T* create()" and "const String getProductName()" (see DefaultParamHandler for details)
  */
  template<UInt D>
  class BaseModel : public DefaultParamHandler
  {
  public:
    typedef double IntensityType;
    typedef double CoordinateType;
    typedef DPosition<D> PositionType;
    typedef typename DPeak<D>::Type PeakType;
    typedef std::vector<PeakType> SamplesType;


    /// Default constructor.
    BaseModel() : DefaultParamHandler("BaseModel")
    {
      defaults_.setValue("cutoff", 0.0, "Low intensity cutoff of the model.  Peaks below this intensity are not considered part of the model.");
    }

    /// copy constructor
    BaseModel(const BaseModel& source) : DefaultParamHandler(source), cut_off_(source.cut_off_)
    {
    }

    /// Destructor
    ~BaseModel() override
    {
    }

    /// assignment operator
    BaseModel& operator=(const BaseModel& source)
    {
      if (&source == this)
        return *this;

      DefaultParamHandler::operator=(source);
      cut_off_ = source.cut_off_;

      return *this;
    }

    /// register all derived classes here (implemented in file BaseModel_impl.h)
    static void registerChildren();

    /// access model predicted intensity at position @p pos
    virtual IntensityType getIntensity(const PositionType& pos) const = 0;

    /// check if position @p pos is part of the model regarding the models cut-off.
    virtual bool isContained(const PositionType& pos) const
    {
      return getIntensity(pos) >= cut_off_;
    }

    /**@brief Convenience function to set the intensity of a peak to the
    predicted intensity at its current position, calling virtual void
    getIntensity().
    */
    template<typename PeakType>
    void fillIntensity(PeakType& peak) const
    {
      peak.setIntensity(getIntensity(peak.getPosition()));
    }

    /**@brief Convenience function that applies fillIntensity() to an iterator
    range.
    */
    template<class PeakIterator>
    void fillIntensities(PeakIterator begin, PeakIterator end) const
    {
      for (PeakIterator it = begin; it != end; ++it)
      {
        fillIntensity(*it);
      }
    }

    /// get cutoff value
    virtual IntensityType getCutOff() const
    {
      return cut_off_;
    }

    /// set cutoff value
    virtual void setCutOff(IntensityType cut_off)
    {
      cut_off_ = cut_off;
      param_.setValue("cutoff", cut_off_);
    }

    /// get reasonable set of samples from the model (i.e. for printing)
    virtual void getSamples(SamplesType& cont) const = 0;

    /// fill stream with reasonable set of samples from the model (i.e. for printing)
    virtual void getSamples(std::ostream& os)
    {
      SamplesType samples;
      getSamples(samples);
      for (typename SamplesType::const_iterator it = samples.begin(); it != samples.end(); ++it)
      {
        os << *it << std::endl;
      }
    }

  protected:
    IntensityType cut_off_;

    void updateMembers_() override
    {
      cut_off_ = (double)param_.getValue("cutoff");
    }
  };
} // namespace OpenMS
