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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_PARENTPEAKMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_PARENTPEAKMOWER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <vector>

namespace OpenMS
{

  /**
    @brief ParentPeakMower gets rid of high peaks that could stem from unfragmented precursor ions

    @htmlinclude OpenMS_ParentPeakMower.parameters

    @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI ParentPeakMower :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    ParentPeakMower();

    /// copy constructor
    ParentPeakMower(const ParentPeakMower& source);

    /// destructor
    ~ParentPeakMower() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    ParentPeakMower& operator=(const ParentPeakMower& source);
    // @}

    // @name Accessors
    // @{
    ///

    ///
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType& spectrum)
    {
      typedef typename SpectrumType::Iterator Iterator;

      clean_all_charge_states_ = (Int)param_.getValue("clean_all_charge_states");
      consider_NH3_loss_ = (Int)param_.getValue("consider_NH3_loss");
      consider_H2O_loss_ = (Int)param_.getValue("consider_H2O_loss");
      window_size_ = (double)param_.getValue("window_size");
      reduce_by_factor_ = (Int)param_.getValue("reduce_by_factor");
      factor_ = (double)param_.getValue("factor");
      set_to_zero_ = (Int)param_.getValue("set_to_zero");

      if (spectrum.getMSLevel() == 1)
      {
        std::cerr << "Error: ParentPeakMower cannot be applied to MS level 1" << std::endl;
        return;
      }

      //get precursor peak position precursor peak
      double pre_pos = 0.0;
      if (!spectrum.getPrecursors().empty()) pre_pos = spectrum.getPrecursors()[0].getMZ();

      if (pre_pos == 0)
      {
        std::cerr << "ParentPeakMower: Warning, Precursor Position not set" << std::endl;
        return;
      }

      Size pre_charge = spectrum.getPrecursors()[0].getCharge();
      if (pre_charge == 0)
      {
        default_charge_ = (Size)param_.getValue("default_charge");
        std::cerr << "ParentPeakMower: Warning, Precursor charge not set, assuming default charge (" << default_charge_ << ")" << std::endl;
        pre_charge = default_charge_;
      }

      pre_pos *= pre_charge;

      // identify the ranges which are to be considered
      std::vector<DRange<1> > ranges;
      for (Size z = 1; z <= pre_charge; ++z)
      {
        if (clean_all_charge_states_ || z == pre_charge)
        {
          // no adjusting needed for this charge
          DPosition<1> pre_z_pos, pos;
          DRange<1> range;

          // adjust the m/z by weight of precursor and charge
          pre_z_pos = DPosition<1>(pre_pos / double(z));
          range = DRange<1>(pre_z_pos - window_size_, pre_z_pos + window_size_);
          ranges.push_back(range);

          if (consider_NH3_loss_)
          {
            pos = DPosition<1>(pre_z_pos - 17.0 / double(z));
            range = DRange<1>(pos - window_size_, pos + window_size_);
            ranges.push_back(range);
          }
          if (consider_H2O_loss_)
          {
            pos = DPosition<1>(pre_z_pos - 18.0 / double(z));
            range = DRange<1>(pos - window_size_, pos + window_size_);
            ranges.push_back(range);
          }
        }
      }

//for (std::vector<DRange<1> >::const_iterator rit = ranges.begin(); rit != ranges.end(); ++rit)
//{
//std::cerr << *rit << std::endl;
//}

      // apply the intensity reduction to the collected ranges
      for (Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
      {
        for (std::vector<DRange<1> >::const_iterator rit = ranges.begin(); rit != ranges.end(); ++rit)
        {
          if (rit->encloses(it->getPosition()))
          {
            if (reduce_by_factor_)
            {
              it->setIntensity(it->getIntensity() / factor_);
              break;
            }

            if (set_to_zero_)
            {
              it->setIntensity(0.0);
              break;
            }
          }
        }
      }

      return;
    }

    void filterPeakSpectrum(PeakSpectrum& spectrum);

    void filterPeakMap(PeakMap& exp);

    //TODO reimplement DefaultParamHandler::updateMembers_()

    //@}

private:
    Size default_charge_;
    bool clean_all_charge_states_;
    bool consider_NH3_loss_;
    bool consider_H2O_loss_;
    double window_size_;
    bool reduce_by_factor_;
    double factor_;
    bool set_to_zero_;

  };

}
#endif // OPENMS_FILTERING/TRANSFORMERS_PARENTPEAKMOWER_H
