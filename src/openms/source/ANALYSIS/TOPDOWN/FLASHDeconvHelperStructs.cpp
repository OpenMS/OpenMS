// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>

namespace OpenMS
{
  FLASHDeconvHelperStructs::PrecalculatedAveragine::PrecalculatedAveragine(const double min_mass,
                                                                           const double max_mass,
                                                                           const double delta,
                                                                           CoarseIsotopePatternGenerator *generator,
                                                                           const bool use_RNA_averagine)
      :
      mass_interval_(delta), min_mass_(min_mass)
  {
    int i = 0;
    while (true)
    {
      double mass_delta = i * mass_interval_;
      i++;
      if (mass_delta < min_mass)
      {
        continue;
      }
      if (mass_delta > max_mass)
      {
        break;
      }
      auto iso = use_RNA_averagine ? generator->estimateFromRNAWeight(mass_delta) : generator->estimateFromPeptideWeight(mass_delta);
      double factor = .01;
      iso.trimRight(factor * iso.getMostAbundant().getIntensity());

      double norm = .0;
      Size most_abundant_index = 0;
      double most_abundant_int = 0;

      for (Size k = 0; k < iso.size(); k++)
      {
        norm += iso[k].getIntensity() * iso[k].getIntensity();
        if (most_abundant_int >= iso[k].getIntensity())
        {
          continue;
        }
        most_abundant_int = iso[k].getIntensity();
        most_abundant_index = k;
      }

      Size left_index = most_abundant_index;
      for (Size k = 0; k <= most_abundant_index; k++)
      {
        if (iso[k].getIntensity() > most_abundant_int * factor)
        {
          break;
        }
        norm -= iso[k].getIntensity() * iso[k].getIntensity();
        left_index--;
        iso[k].setIntensity(0);
      }

      Size right_index = iso.size() - 1 - most_abundant_index;
      /*for (Size k = iso.size() - 1; k >= most_abundant_index; k--)
      {
        if (iso[k].getIntensity() > most_abundant_int * factor)
        {
          break;
        }
        norm -= iso[k].getIntensity() * iso[k].getIntensity();
        right_index--;
        iso[k].setIntensity(0);
      }*/

      isotope_end_indices_.push_back(right_index);
      isotope_start_indices_.push_back(left_index);
      average_mono_mass_difference_.push_back(iso.averageMass() - iso[0].getMZ());
      norms_.push_back(norm);
      isotopes_.push_back(iso);
    }
  }

  IsotopeDistribution FLASHDeconvHelperStructs::PrecalculatedAveragine::get(const double mass) const
  {
    Size i = (Size) (.5 + std::max(.0, mass - min_mass_) / mass_interval_);
    i = i >= isotopes_.size() ? isotopes_.size() - 1 : i;
    return isotopes_[i];
  }

  int FLASHDeconvHelperStructs::PrecalculatedAveragine::getMaxIsotopeIndex() const
  {
    return max_isotope_index_;
  }


  double FLASHDeconvHelperStructs::PrecalculatedAveragine::getNorm(const double mass) const
  {
    Size i = (Size) (.5 + std::max(.0, mass - min_mass_) / mass_interval_);
    i = i >= isotopes_.size() ? isotopes_.size() - 1 : i;
    return norms_[i];
  }

  Size FLASHDeconvHelperStructs::PrecalculatedAveragine::getIsotopeStartIndex(const double mass) const
  {
    Size i = (Size) (.5 + std::max(.0, mass - min_mass_) / mass_interval_);
    i = i >= isotopes_.size() ? isotopes_.size() - 1 : i;
    return isotope_start_indices_[i];
  }

  double FLASHDeconvHelperStructs::PrecalculatedAveragine::getAverageMassDelta(const double mass) const
  {
    Size i = (Size) (.5 + std::max(.0, mass - min_mass_) / mass_interval_);
    i = i >= isotopes_.size() ? isotopes_.size() - 1 : i;
    return average_mono_mass_difference_[i];
  }

  Size FLASHDeconvHelperStructs::PrecalculatedAveragine::getIsotopeEndIndex(const double mass) const
  {
    Size i = (Size) (.5 + std::max(.0, mass - min_mass_) / mass_interval_);
    i = i >= isotopes_.size() ? isotopes_.size() - 1 : i;
    return isotope_end_indices_[i];
  }

  void FLASHDeconvHelperStructs::PrecalculatedAveragine::setMaxIsotopeIndex(const int index)
  {
    max_isotope_index_ = index;
  }

  FLASHDeconvHelperStructs::LogMzPeak::LogMzPeak(const Peak1D& peak, const bool positive) :
      mz(peak.getMZ()),
      intensity(peak.getIntensity()),
      logMz(getLogMz(peak.getMZ(), positive)),
      abs_charge(0),
      is_positive(positive),
      isotopeIndex(0)
  {
  }

  double FLASHDeconvHelperStructs::LogMzPeak::getUnchargedMass()
  {
    if (abs_charge == 0)
    {
      return .0;
    }
    if (mass <= 0)
    {
      mass = (mz - getChargeMass(is_positive)) * abs_charge;
    }
    return mass;
  }

  bool FLASHDeconvHelperStructs::LogMzPeak::operator<(const LogMzPeak& a) const
  {
    return this->logMz < a.logMz;
  }

  bool FLASHDeconvHelperStructs::LogMzPeak::operator>(const LogMzPeak& a) const
  {
    return this->logMz > a.logMz;
  }

  bool FLASHDeconvHelperStructs::LogMzPeak::operator==(const LogMzPeak& a) const
  {
    return this->logMz == a.logMz;
  }

  class LogMzPeakHashFunction
  {
  public:

    size_t operator()(const FLASHDeconvHelperStructs::LogMzPeak& peak) const
    {
      return std::hash<double>()(peak.mz);
    }
  };

  FLASHDeconvHelperStructs::PrecalculatedAveragine FLASHDeconvHelperStructs::calculateAveragines(const double max_mass,
                                                                                                 const bool use_RNA_averagine)
  {
    auto generator = new CoarseIsotopePatternGenerator();

    auto iso = use_RNA_averagine ?
               generator->estimateFromRNAWeight(max_mass) :
               generator->estimateFromPeptideWeight(max_mass);
    iso.trimRight(0.01 * iso.getMostAbundant().getIntensity());

    generator->setMaxIsotope(iso.size());
    auto avg = FLASHDeconvHelperStructs::PrecalculatedAveragine(50, max_mass, 25, generator, use_RNA_averagine);
    avg.setMaxIsotopeIndex(iso.size() - 1);
    return avg;
  }

  double FLASHDeconvHelperStructs::getChargeMass(const bool positive)
  {
    return (positive ? Constants::PROTON_MASS_U : Constants::ELECTRON_MASS_U);
  }


  double FLASHDeconvHelperStructs::getLogMz(const double mz, const bool positive)
  {
    return std::log(mz - getChargeMass(positive));
  }
}
