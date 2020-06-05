// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>

using namespace std;

namespace OpenMS
{

  MultiplexIsotopicPeakPattern::MultiplexIsotopicPeakPattern(int c, int ppp, MultiplexDeltaMasses ms, int msi) :
    charge_(c), peaks_per_peptide_(ppp), mass_shifts_(ms), mass_shift_index_(msi)
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
