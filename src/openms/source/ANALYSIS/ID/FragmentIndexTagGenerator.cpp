// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGenerator.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/QC/QCBase.h>
#include <functional>


using namespace std;


namespace OpenMS
{
  std::vector<std::string> aminoAcids = {
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
  };
  std::vector<float> monoisotopic_masses = {
    71.03711, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891, 113.08406, 128.09496, 131.04049, 114.04293, 97.05276, 128.05858, 156.10111, 87.03203, 101.04768, 99.06841, 186.07932, 163.06332
  };


  FragmentIndexTagGenerator::FragmentIndexTagGenerator(const OpenMS::MSSpectrum& spectrum)
  {
    spectrum_ = spectrum;
    selected_peaks_.resize(spectrum_.size());

  }
  FragmentIndexTagGenerator::~FragmentIndexTagGenerator() = default;

  FragmentIndexTagGenerator::FragmentIndexTagGenerator(const OpenMS::FragmentIndexTagGenerator& cp) = default;

  FragmentIndexTagGenerator& FragmentIndexTagGenerator::operator=(const OpenMS::FragmentIndexTagGenerator& source)
  {
    if(this != &source)
    {
      spectrum_ = source.spectrum_;
      selected_peaks_ = source.selected_peaks_;
      n = source.n;
    }
    return *this;
  }


  void FragmentIndexTagGenerator::setMSSpectrum(const MSSpectrum& spectrum)
  {
    spectrum_.clear(true);
    selected_peaks_.clear();

    spectrum_ = spectrum;
    selected_peaks_.resize(spectrum_.size());
    n = spectrum_.getPrecursors()[0].getMZ()/ 35;  // TODO: Does this make sense??
  }


  void FragmentIndexTagGenerator::globalSelection()
  {
    vector<IdxAndIntensity> idx_intensity;
    uint32_t idx_count = 0;
    for(Peak1D peak: spectrum_){
      idx_intensity.emplace_back(idx_count, peak.getIntensity());
      idx_count++;
    }

    std::partial_sort(idx_intensity.begin(), idx_intensity.begin() + n, idx_intensity.end(), [](const IdxAndIntensity a, const IdxAndIntensity b)
                      {
                        if(a.intensity_ == b.intensity_)
                          return a.idx_ < b.idx_;
                        else
                          return a.idx_ > b.idx_;
                      }

    );

    for(auto iai = idx_intensity.begin();iai != idx_intensity.end() && iai != idx_intensity.begin() + n; iai++){
      selected_peaks_[iai->idx_] = true;

    }

  }

  void FragmentIndexTagGenerator::localSelection()
  {
    size_t lower_peak = 0;
    size_t upper_peak = 0;
    for(float window_da = 0; window_da < spectrum_.getPrecursors()[0].getUnchargedMass(); window_da += 35.0){ //window increment = 35 Da
      while (upper_peak < spectrum_.size() && spectrum_[upper_peak].getMZ() <= (window_da +70)) upper_peak++;   //window size      = 70 Da
      size_t sum = 0; // how many peaks are in the 70-Da window
      vector<IdxAndIntensity> temp_snippet; // collect all idx within the window and sort them afterward
      for(auto sele_iter = selected_peaks_.begin() + lower_peak; sele_iter != selected_peaks_.begin()+upper_peak; sele_iter++){
        sum += static_cast<size_t>(*sele_iter);
        size_t current_peak =sele_iter - selected_peaks_.begin();
        temp_snippet.emplace_back(current_peak, spectrum_[current_peak].getIntensity());
      }
      sort(temp_snippet.begin(), temp_snippet.end(),[](IdxAndIntensity a, IdxAndIntensity b){return (a.intensity_ > b.intensity_);});
      while((sum < 2) && (sum < temp_snippet.size())){                 // now select new peaks if needed for this window
        selected_peaks_[temp_snippet[sum].idx_] = true;
        sum++;
      }

      lower_peak = upper_peak;
    }
  }


  void FragmentIndexTagGenerator::generateAllMultiPeaksFast(std::vector<MultiPeak>& all_multi_peaks,
                                                            size_t depth,
                                                            float tolerance,
                                                            bool tolerance_unit)
  {
    //first check if spectrum is sorted!!!
    if (!spectrum_.isSorted())
    {
      OPENMS_LOG_WARN << "Spectrum was not sorted. Will be sorted now!" << endl;
      spectrum_.sortByPosition();
    }
    for (size_t i = 0; i < spectrum_.size(); i++)  //loop over all peaks. and execute the recursion
    {
      // initialize the MultiPeak
      float actual_tol;
      if(tolerance_unit)
        actual_tol =  Math::ppmToMass(tolerance, (float) spectrum_[i].getMZ());
      else
        actual_tol = tolerance;
      MultiPeak multi_peak((float) spectrum_[i].getMZ(),(float) spectrum_[i].getIntensity());
      generateAllMultiPeaksFastRecursion(all_multi_peaks, multi_peak, i, depth, actual_tol);

    }
  }

  void FragmentIndexTagGenerator::generateAllMultiPeaksFastRecursion(std::vector<MultiPeak>& all_multi_peaks, FragmentIndexTagGenerator::MultiPeak& multi_peak, size_t current_peak_idx,
                                                                     size_t recursion_step, float tolerance)
  {
      if (recursion_step == 0)
      {
        all_multi_peaks.push_back(multi_peak);
        return;
      }
      size_t j = current_peak_idx + 1;
      while (j < spectrum_.size())
      {
        auto delta_mz = static_cast<float>(spectrum_[j].getMZ() - spectrum_[current_peak_idx].getMZ());
        if (delta_mz > 190)
        {
          break;
        }
        for (size_t mim = 0; mim < monoisotopic_masses.size(); mim++)
        {
          if ((monoisotopic_masses[mim] - tolerance <= delta_mz) &&(delta_mz <= monoisotopic_masses[mim] + tolerance)){
            FragmentIndexTagGenerator::MultiPeak temp_copy(multi_peak);
            temp_copy.addFollowUpPeak(delta_mz, aminoAcids[mim]);
            temp_copy.addScore(spectrum_[j].getIntensity());
            generateAllMultiPeaksFastRecursion(all_multi_peaks, temp_copy, j, recursion_step - 1, tolerance);
            break;
          }
        }
        j++;
      }
  }


  void FragmentIndexTagGenerator::generateAllMultiFragments(std::vector<FragmentIndexTagGenerator::MultiFragment>& multi_frags,
                                               size_t depth,
                                               size_t peptide_idx, float frag_min_mz, float frag_max_mz)
  {
    if(spectrum_.getStringDataArrays().empty())
    {
      OPENMS_LOG_WARN << "The provided spectrum has no ion-type info" << endl;
      return;
    }
    const PeakSpectrum::StringDataArray  & ion_types = spectrum_.getStringDataArrays().at(0);
    for(size_t i = 0; i < spectrum_.size(); i++){                   // loop through all peaks, check if they have mz in the window
      if((frag_min_mz > spectrum_[i].getMZ()) && (spectrum_[i].getMZ() > frag_max_mz))
        continue;
      vector<float> temp_follow_up;
      size_t j = i +1;
      size_t last_j = i;
      while(temp_follow_up.size()< depth && j < spectrum_.size()){               // look at all follow up peaks of the same ion type
        if(ion_types[i].substr(0,1) == ion_types[j].substr(0,1)){
          temp_follow_up.push_back((float)(spectrum_[j].getMZ() - spectrum_[last_j].getMZ()));
          last_j = j;
        }
        j++;
      }
      if(temp_follow_up.size() == depth)
        multi_frags.emplace_back(peptide_idx, (float)spectrum_[i].getMZ(), temp_follow_up);
    }
  }


  FragmentIndexTagGenerator::MultiFragment::MultiFragment() :
      peptide_idx_(0),
      fragment_mz_(0),
      follow_up_peaks_(3,0.0)
  {
  }

  FragmentIndexTagGenerator::MultiFragment::MultiFragment(OpenMS::UInt32 peptide_idx, float fragment_mz, const vector<float>& follow_up):
      peptide_idx_(peptide_idx),
      fragment_mz_(fragment_mz),
      follow_up_peaks_(follow_up)
  //follow_up_peaks_AA_(std::move(sequ))
  {
  }

  FragmentIndexTagGenerator::MultiFragment::MultiFragment(OpenMS::UInt32 peptide_idx, float fragment_mz, const FragmentIndexTagGenerator::MultiPeak& multiPeak) :
      peptide_idx_(peptide_idx), fragment_mz_(fragment_mz)
  {
    follow_up_peaks_ = multiPeak.getFollowUpPeaks();
    //follow_up_peaks_AA_ = multiPeak.getFollowUpPeaksAa();
  }

  FragmentIndexTagGenerator::MultiFragment::MultiFragment(const FragmentIndexTagGenerator::MultiFragment& other)  = default;

  FragmentIndexTagGenerator::MultiFragment& FragmentIndexTagGenerator::MultiFragment::operator=(const FragmentIndexTagGenerator::MultiFragment& other)
  {
    if(&other == this)
      return *this;
    peptide_idx_ = other.peptide_idx_;
    fragment_mz_ = other.fragment_mz_;
    follow_up_peaks_ = other.follow_up_peaks_;
    //follow_up_peaks_AA_ = other.follow_up_peaks_AA_;
    return *this;
  }

  void FragmentIndexTagGenerator::MultiFragment::swap(FragmentIndexTagGenerator::MultiFragment& other)
  {
    follow_up_peaks_.swap(other.follow_up_peaks_);
    //follow_up_peaks_AA_.swap(other.follow_up_peaks_AA_);
    std::swap(peptide_idx_, other.peptide_idx_);
    std::swap(fragment_mz_, other.fragment_mz_);
  }

  UInt32 FragmentIndexTagGenerator::MultiFragment::getPeptideIdx() const
  {
    return peptide_idx_;
  }
  float FragmentIndexTagGenerator::MultiFragment::getFragmentMz() const
  {
    return fragment_mz_;
  }

  const vector<float>& FragmentIndexTagGenerator::MultiFragment::getFollowUpPeaks() const
  {
    return follow_up_peaks_;
  }

  FragmentIndexTagGenerator::MultiPeak::MultiPeak() : peak_(), score_(0), follow_up_peaks({})
  {
  }
  FragmentIndexTagGenerator::MultiPeak::MultiPeak(float peak, float score) : peak_(peak), score_(score), follow_up_peaks()
  {
  }
  FragmentIndexTagGenerator::MultiPeak::MultiPeak(const FragmentIndexTagGenerator::MultiPeak& other)    :
      peak_(other.peak_),
      score_(other.score_),
      follow_up_peaks_AA(other.follow_up_peaks_AA),
      follow_up_peaks(other.follow_up_peaks) {
      };
  FragmentIndexTagGenerator::MultiPeak& FragmentIndexTagGenerator::MultiPeak::operator=(const FragmentIndexTagGenerator::MultiPeak& other)
  {
    if (&other == this)
      return *this;
    peak_ = other.peak_;
    score_ = other.score_;
    follow_up_peaks_AA = other.follow_up_peaks_AA;
    follow_up_peaks = other.follow_up_peaks;
    return  *this;
  }
  void FragmentIndexTagGenerator::MultiPeak::addScore(float score)
  {
    score_ += score;
  }

  float FragmentIndexTagGenerator::MultiPeak::getPeak() const
  {
    return peak_;
  }
  float FragmentIndexTagGenerator::MultiPeak::getScore() const
  {
    return score_;
  }
  const string& FragmentIndexTagGenerator::MultiPeak::getFollowUpPeaksAa() const
  {
    return follow_up_peaks_AA;
  }
  const vector<float>& FragmentIndexTagGenerator::MultiPeak::getFollowUpPeaks() const
  {
    return follow_up_peaks;
  }

  void FragmentIndexTagGenerator::MultiPeak::addFollowUpPeak(float distance, const std::string &AA)
  {
    follow_up_peaks_AA += AA;
    follow_up_peaks.push_back(distance);
  }

}
