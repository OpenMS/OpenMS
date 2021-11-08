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
#include <sstream>

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
      double mass = i * mass_interval_;
      i++;
      if (mass < min_mass)
      {
        continue;
      }
      if (mass > max_mass)
      {
        break;
      }
      auto iso = use_RNA_averagine ?
                 generator->estimateFromRNAWeight(mass) :
                 generator->estimateFromPeptideWeight(mass);

      const double min_pwr = .999;
      const Size min_iso_length = 3;
      double total_pwr = .0;
      int most_abundant_index_ = 0;
      double most_abundant_int = 0;

      for (Size k = 0; k < iso.size(); k++)
      {
        total_pwr += iso[k].getIntensity() * iso[k].getIntensity();
        if (most_abundant_int >= iso[k].getIntensity())
        {
          continue;
        }
        most_abundant_int = iso[k].getIntensity();
        most_abundant_index_ = k;
      }

      int left_count = 0;
      int right_count = iso.size() - 1;
      int trim_count = 0;
      while (iso.size() - trim_count > min_iso_length)
      {
        double lint = iso[left_count].getIntensity();
        double rint = iso[right_count].getIntensity();
        double pwr;
        bool trim_left = true;
        if (lint < rint)
        {
          pwr = lint * lint;
        }
        else
        {
          pwr = rint * rint;
          trim_left = false;
        }
        if (total_pwr - pwr < total_pwr * min_pwr)
        {
          break;
        }
        total_pwr -= pwr;
        trim_count++;
        if (trim_left)
        {
          iso[left_count].setIntensity(0);
          left_count++;
        }
        else
        {
          iso[right_count].setIntensity(0); // for trimming
          right_count--;
        }
      }
      left_count = most_abundant_index_ - left_count;
      right_count = right_count - most_abundant_index_;
      iso.trimRight(1e-10);

      for (Size k = 0; k < iso.size(); k++)
      {
        double ori_int = iso[k].getIntensity();
        iso[k].setIntensity(ori_int / sqrt(total_pwr));
      }
      left_count = left_count < 0 ? 0 : left_count;
      right_count = right_count < 0 ? 0 : right_count;

      apex_index_.push_back(most_abundant_index_);
      right_count_from_apex_.push_back(1 + right_count);
      left_count_from_apex_.push_back(1 + left_count);
      average_mono_mass_difference_.push_back(iso.averageMass() - iso[0].getMZ());
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


  //double FLASHDeconvHelperStructs::PrecalculatedAveragine::getNorm(const double mass) const {
  //    Size i = (Size) (.5 + std::max(.0, mass - min_mass_) / mass_interval_);
  //    i = i >= isotopes_.size() ? isotopes_.size() - 1 : i;
  //    return norms_[i];
  //}

  Size FLASHDeconvHelperStructs::PrecalculatedAveragine::getLeftCountFromApex(const double mass) const
  {
    Size i = (Size) (.5 + std::max(.0, mass - min_mass_) / mass_interval_);
    i = i >= isotopes_.size() ? isotopes_.size() - 1 : i;
    return (Size) left_count_from_apex_[i];
  }

  double FLASHDeconvHelperStructs::PrecalculatedAveragine::getAverageMassDelta(const double mass) const
  {
    Size i = (Size) (.5 + std::max(.0, mass - min_mass_) / mass_interval_);
    i = i >= isotopes_.size() ? isotopes_.size() - 1 : i;
    return average_mono_mass_difference_[i];
  }

  Size FLASHDeconvHelperStructs::PrecalculatedAveragine::getRightCountFromApex(const double mass) const
  {
    Size i = (Size) (.5 + std::max(.0, mass - min_mass_) / mass_interval_);
    i = i >= isotopes_.size() ? isotopes_.size() - 1 : i;
    return (Size) right_count_from_apex_[i];
  }

  Size FLASHDeconvHelperStructs::PrecalculatedAveragine::getApexIndex(const double mass) const
  {
    Size i = (Size) (.5 + std::max(.0, mass - min_mass_) / mass_interval_);
    i = i >= isotopes_.size() ? isotopes_.size() - 1 : i;
    return apex_index_[i];
  }

  void FLASHDeconvHelperStructs::PrecalculatedAveragine::setMaxIsotopeIndex(const int index)
  {
    max_isotope_index_ = index;
  }

  FLASHDeconvHelperStructs::LogMzPeak::LogMzPeak(const Peak1D &peak, const bool positive) :
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

  bool FLASHDeconvHelperStructs::LogMzPeak::operator<(const LogMzPeak &a) const
  {
    return this->logMz < a.logMz;
  }

  bool FLASHDeconvHelperStructs::LogMzPeak::operator>(const LogMzPeak &a) const
  {
    return this->logMz > a.logMz;
  }

  bool FLASHDeconvHelperStructs::LogMzPeak::operator==(const LogMzPeak &a) const
  {
    return this->logMz == a.logMz;
  }

  FLASHDeconvHelperStructs::PrecalculatedAveragine
  FLASHDeconvHelperStructs::calculateAveragines(const double max_mass,
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
    return (positive ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U);
  }


  double FLASHDeconvHelperStructs::getLogMz(const double mz, const bool positive)
  {
    return std::log(mz - getChargeMass(positive));
  }

  FLASHDeconvHelperStructs::TopPicItem::TopPicItem(String in)
  {
    str_ = in;
    std::vector<String> results;
    std::stringstream tmp_stream(in);
    String str;
    //Data file name	Prsm ID	Spectrum ID	Fragmentation	Scan(s)	Retention time	#peaks	Charge	Precursor mass
    // Adjusted precursor mass	Proteoform ID	Feature intensity	Feature score	Protein accession	Protein description
    // First residue	Last residue	Proteoform	#unexpected modifications	MIScore	#variable PTMs	#matched peaks
    // #matched fragment ions	E-value	Spectrum-level Q-value	Proteoform-level Q-value

    while (getline(tmp_stream, str, '\t'))
    {
      results.push_back(str);
    }
    prsm_id_ = stoi(results[1]);
    spec_id_ = stoi(results[2]);
    scan_ = stoi(results[4]);
    rt_ = stod(results[5]);
    peak_count_ = stoi(results[6]);
    charge_ = stoi(results[7]);
    precursor_mass_ = stod(results[8]);
    adj_precursor_mass_ = stod(results[9]);
    proteform_id_ = stoi(results[10]);
    if (results[11].hasPrefix("-"))
    {
      intensity_ = 0;
    }
    else
    {
      intensity_ = stod(results[11]);
    }
    String acc = results[13];
    int first = acc.find("|");
    int second = acc.find("|", first + 1);
    protein_acc_ = acc.substr(first + 1, second - first - 1);
    first_residue_ = stoi(results[15]);
    last_residue_ = stoi(results[16]);
    if (stoi(results[18]) == 0)
    {
      //unexp_mod_ = .0;
    }
    else
    {
      String seq = results[17];
      int loc = 0;
      //int off = seq.find(".", 0);
      while (seq.find("[", loc) != String::npos)
      {
        // mod_first_.push_back(seq.find("(", loc) - off -1);
        // mod_last_.push_back(seq.find(")", loc) - off -3);
        loc = seq.find("[", loc);
        String sub = seq.substr(loc + 1, seq.find("]", loc) - 1 - loc);
        double mmass = .0;
        if (isdigit(sub[sub.length() - 1]))
        {
          mmass = stod(sub);
        }
        else if (sub == "Acetyl")
        {
          mmass = 42.010565;
        }
        else if (sub == "Phospho")
        {
          mmass = 79.966331;
        }
        else if (sub == "Oxidation")
        {
          mmass = 15.994915;
        }
        else if (sub == "Methyl")
        {
          mmass = 14.015650;
        }
        unexp_mod_.push_back(mmass);
        loc++;
      }
    }


    matched_peaks_ = stoi(results[21]);
    matched_frags_ = stoi(results[22]);
    e_value_ = stod(results[23]);
    if (results[24] != "-")
    {
      spec_q_value_ = stod(results[24]);
      proteofrom_q_value_ = stod(results[25]);
    }
  }

  bool FLASHDeconvHelperStructs::TopPicItem::operator<(const FLASHDeconvHelperStructs::TopPicItem &a) const
  {
    return this->scan_ < a.scan_;
  }

  bool FLASHDeconvHelperStructs::TopPicItem::operator>(const FLASHDeconvHelperStructs::TopPicItem &a) const
  {
    return this->scan_ > a.scan_;
  }

  bool FLASHDeconvHelperStructs::TopPicItem::operator==(const FLASHDeconvHelperStructs::TopPicItem &other) const
  {
    return this->scan_ == other.scan_;
  }
}
