// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/DPeak.h>

namespace OpenMS
{

  /**
    @brief Abstract base class for all D-dimensional models.

    Every derived class has to implement the static functions
    "T* create()" and "const String getProductName()" (see DefaultParamHandler for details)
  */
  template <UInt D>
  class BaseModel :
    public DefaultParamHandler
  {

public:

    typedef double IntensityType;
    typedef double CoordinateType;
    typedef DPosition<D> PositionType;
    typedef typename DPeak<D>::Type PeakType;
    typedef std::vector<PeakType> SamplesType;


    /// Default constructor.
    BaseModel() :
      DefaultParamHandler("BaseModel")
    {
      defaults_.setValue("cutoff", 0.0, "Low intensity cutoff of the model.  Peaks below this intensity are not considered part of the model.");
    }

    /// copy constructor
    BaseModel(const BaseModel & source) :
      DefaultParamHandler(source),
      cut_off_(source.cut_off_)
    {
    }

    /// Destructor
    ~BaseModel() override
    {
    }

    /// assignment operator
    virtual BaseModel & operator=(const BaseModel & source)
    {
      if (&source == this) return *this;

      DefaultParamHandler::operator=(source);
      cut_off_ = source.cut_off_;

      return *this;
    }

    /// register all derived classes here (implemented in file BaseModel_impl.h)
    static void registerChildren();

    /// access model predicted intensity at position @p pos
    virtual IntensityType getIntensity(const PositionType & pos) const = 0;

    /// check if position @p pos is part of the model regarding the models cut-off.
    virtual bool isContained(const PositionType & pos) const
    {
      return getIntensity(pos) >= cut_off_;
    }

    /**@brief Convenience function to set the intensity of a peak to the
    predicted intensity at its current position, calling virtual void
    getIntensity().
    */
    template <typename PeakType>
    void fillIntensity(PeakType & peak) const
    {
      peak.setIntensity(getIntensity(peak.getPosition()));
    }

    /**@brief Convenience function that applies fillIntensity() to an iterator
    range.
    */
    template <class PeakIterator>
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
    virtual void getSamples(SamplesType & cont) const = 0;

    /// fill stream with reasonable set of samples from the model (i.e. for printing)
    virtual void getSamples(std::ostream & os)
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
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H
