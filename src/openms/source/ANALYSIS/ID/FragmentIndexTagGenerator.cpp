// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGenerator.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGeneratorNode.h>
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
  TagGenerator::TagGenerator(const OpenMS::MSSpectrum& spectrum)
  {
    spectrum_ = spectrum;
    selected_peaks_.resize(spectrum_.size());
    n = spectrum_.getPrecursors()[0].getMZ()/ 35;  // TODO: Does this make sense??
  }
  TagGenerator::~TagGenerator() = default;

  TagGenerator::TagGenerator(const OpenMS::TagGenerator& cp) = default;

  TagGenerator& TagGenerator::operator=(const OpenMS::TagGenerator& source)
  {
    if(this != &source)
    {
      spectrum_ = source.spectrum_;
      selected_peaks_ = source.selected_peaks_;
      n = source.n;
      dag_ = source.dag_;
    }
    return *this;
  }


  void TagGenerator::setMSSpectrum(const MSSpectrum& spectrum)
  {
    spectrum_.clear(true);
    selected_peaks_.clear();
    dag_.clear();

    spectrum_ = spectrum;
    selected_peaks_.resize(spectrum_.size());
    n = spectrum_.getPrecursors()[0].getMZ()/ 35;  // TODO: Does this make sense??
  }


  void TagGenerator::globalSelection()
  {
    vector<IdxAndIntensity> idx_intensity;
    uint32_t idx_count = 0;
    for(Peak1D peak: spectrum_){
      idx_intensity.emplace_back(idx_count, peak.getIntensity());
      idx_count++;
    }
    //we can use the same heap structure as in the scoring to get the top k hits
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

  void TagGenerator::localSelection()
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

  void TagGenerator::generateAllNodes(uint32_t max_charge)
  {
    //should there already be a graph delete it
    dag_.clear();
    size_t peak_counter = 0;
    for(Peak1D peak: spectrum_){
      for(uint32_t charge = 1; charge <= max_charge; charge++){
        dag_.push_back(make_shared<TagGeneratorNode>(peak, charge, selected_peaks_[peak_counter]));
      }
      peak_counter++;
    }
    sort(dag_.begin(), dag_.end(), [](const shared_ptr<TagGeneratorNode>& a, const shared_ptr<TagGeneratorNode>& b){
      return a->calculateMass() < b->calculateMass();
    });
  }

  void TagGenerator::generateDirectedAcyclicGraph(float fragment_tolerance)
  {
    //first extract the max intensity to later normalize the intensities
    float max = 0;
    for(Peak1D p: spectrum_){
      if(p.getIntensity() > max) max = p.getIntensity();
    }

    for(size_t i = selected_peaks_.size(); i > 0; i--){
      if(selected_peaks_[i-1]){
        shared_ptr<TagGeneratorNode> newNode = make_shared<TagGeneratorNode>(spectrum_[i-1]);
        newNode->calculateConfidence(spectrum_, max);
        for(auto dag_iter = dag_.end() - 1; dag_iter != dag_.begin()-1 ; dag_iter--){    // loop through the complete loop
          if(!newNode->generateConnection(*dag_iter, fragment_tolerance))                      // this function returns false if we out of any
            break;                                                                                  // AA range, in which case we will break the loop
        }
        dag_.push_back(newNode);

      }
    }
  }

  void TagGenerator::generateAllMultiPeaks(std::vector<TagGenerator::MultiPeak>& quad_peaks, size_t depth)
  {
    for(shared_ptr<TagGeneratorNode> node: dag_){
      vector<MultiPeak> quad_peaks_per_node;
      node->generateAllMultiPeaks(quad_peaks_per_node, depth);  // for every node in the dag start a recursiv chain
      quad_peaks.insert(quad_peaks.end(), quad_peaks_per_node.begin(), quad_peaks_per_node.end());
    }
  }

  void TagGenerator::generateAllMultiFragments(std::vector<MultiFragment>& multi_frags,
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


  TagGenerator::MultiFragment::MultiFragment() :
      peptide_idx_(0),
      fragment_mz_(0),
      follow_up_peaks_(3,0.0)
  {
  }

  TagGenerator::MultiFragment::MultiFragment(OpenMS::UInt32 peptide_idx, float fragment_mz, const vector<float>& follow_up):
      peptide_idx_(peptide_idx),
      fragment_mz_(fragment_mz),
      follow_up_peaks_(follow_up)
  //follow_up_peaks_AA_(std::move(sequ))
  {
  }

  TagGenerator::MultiFragment::MultiFragment(OpenMS::UInt32 peptide_idx, float fragment_mz, const OpenMS::MultiPeak& multiPeak) :
      peptide_idx_(peptide_idx), fragment_mz_(fragment_mz)
  {
    follow_up_peaks_ = multiPeak.getFollowUpPeaks();
    //follow_up_peaks_AA_ = multiPeak.getFollowUpPeaksAa();
  }

  TagGenerator::MultiFragment::MultiFragment(const TagGenerator::MultiFragment& other) = default;

  TagGenerator::MultiFragment& TagGenerator::MultiFragment::operator=(const TagGenerator::MultiFragment& other)
  {
    if(&other == this)
      return *this;
    peptide_idx_ = other.peptide_idx_;
    fragment_mz_ = other.fragment_mz_;
    follow_up_peaks_ = other.follow_up_peaks_;
    //follow_up_peaks_AA_ = other.follow_up_peaks_AA_;
    return *this;
  }

  void TagGenerator::MultiFragment::swap(TagGenerator::MultiFragment& other)
  {
    follow_up_peaks_.swap(other.follow_up_peaks_);
    //follow_up_peaks_AA_.swap(other.follow_up_peaks_AA_);
    std::swap(peptide_idx_, other.peptide_idx_);
    std::swap(fragment_mz_, other.fragment_mz_);
  }

  UInt32 TagGenerator::MultiFragment::getPeptideIdx() const
  {
    return peptide_idx_;
  }
  float TagGenerator::MultiFragment::getFragmentMz() const
  {
    return fragment_mz_;
  }

  const vector<float>& TagGenerator::MultiFragment::getFollowUpPeaks() const
  {
    return follow_up_peaks_;
  }

  MultiPeak::MultiPeak() : peak_(), score_(0), follow_up_peaks({})
  {
  }
  MultiPeak::MultiPeak(OpenMS::Peak1D peak, float score) : peak_(peak), score_(score), follow_up_peaks()
  {
  }
  MultiPeak::MultiPeak(const OpenMS::MultiPeak& other) = default;
  MultiPeak& MultiPeak::operator=(const OpenMS::MultiPeak& other)
  {
    if (&other == this)
      return *this;
    peak_ = other.peak_;
    score_ = other.score_;
    follow_up_peaks_AA = other.follow_up_peaks_AA;
    follow_up_peaks = other.follow_up_peaks;
    return  *this;
  }
  void MultiPeak::addScore(float score)
  {
    score_ += score;
  }

  const Peak1D& MultiPeak::getPeak() const
  {
    return peak_;
  }
  float MultiPeak::getScore() const
  {
    return score_;
  }
  const string& MultiPeak::getFollowUpPeaksAa() const
  {
    return follow_up_peaks_AA;
  }
  const vector<float>& MultiPeak::getFollowUpPeaks() const
  {
    return follow_up_peaks;
  }

  void MultiPeak::addFollowUpPeak(float distance, const std::string &AA)
  {
    follow_up_peaks_AA += AA;
    follow_up_peaks.push_back(distance);
  }

}
