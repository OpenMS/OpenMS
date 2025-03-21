// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>

#include <utility>

using namespace std;

namespace OpenMS
{

  MultiplexIsotopicPeakPattern::MultiplexIsotopicPeakPattern(int c, int ppp, MultiplexDeltaMasses ms, int msi) :
    charge_(c), peaks_per_peptide_(ppp), mass_shifts_(std::move(ms)), mass_shift_index_(msi)
  {
    // generate m/z shifts
    for (unsigned i = 0; i < mass_shifts_.getDeltaMasses().size(); ++i)
    {
      for (int j = 0; j < peaks_per_peptide_; ++j)
      {
        const std::vector<MultiplexDeltaMasses::DeltaMass>& delta_masses = mass_shifts_.getDeltaMasses();
        mz_shifts_.push_back((delta_masses[i].delta_mass + j * Constants::C13C12_MASSDIFF_U) / charge_);
      }
    }
  }

  int MultiplexIsotopicPeakPattern::getCharge() const
  {
    return charge_;
  }

  int MultiplexIsotopicPeakPattern::getPeaksPerPeptide() const
  {
    return peaks_per_peptide_;
  }

  MultiplexDeltaMasses MultiplexIsotopicPeakPattern::getMassShifts() const
  {
    return mass_shifts_;
  }

  int MultiplexIsotopicPeakPattern::getMassShiftIndex() const
  {
    return mass_shift_index_;
  }

  unsigned MultiplexIsotopicPeakPattern::getMassShiftCount() const
  {
    return mass_shifts_.getDeltaMasses().size();
  }

  double MultiplexIsotopicPeakPattern::getMassShiftAt(size_t i) const
  {
    return mass_shifts_.getDeltaMasses()[i].delta_mass;
  }

  double MultiplexIsotopicPeakPattern::getMZShiftAt(size_t i) const
  {
    return mz_shifts_[i];
  }

  unsigned MultiplexIsotopicPeakPattern::getMZShiftCount() const
  {
    return mz_shifts_.size();
  }

}
