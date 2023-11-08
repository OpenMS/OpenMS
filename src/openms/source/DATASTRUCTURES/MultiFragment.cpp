//
// Created by trapho on 11/3/23.
//
#include <OpenMS/DATASTRUCTURES/MultiFragment.h>
#include <OpenMS/DATASTRUCTURES/MultiPeak.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/ANALYSIS/ID/TagGeneratorNode.h>

#include <utility>
#include <vector>
#include <functional>

using namespace std;

namespace OpenMS
{


  MultiFragment::MultiFragment() :
  peptide_idx_(0),
  fragment_mz_(0),
  follow_up_peaks_(3,0.0)
  {
  }

  MultiFragment::MultiFragment(OpenMS::Size peptide_idx, double fragment_mz, const vector<double>& follow_up):
  peptide_idx_(peptide_idx),
  fragment_mz_(fragment_mz),
  follow_up_peaks_(follow_up)
  //follow_up_peaks_AA_(std::move(sequ))
  {
  }

  MultiFragment::MultiFragment(OpenMS::Size peptide_idx, double fragment_mz, const OpenMS::MultiPeak& multiPeak) :
  peptide_idx_(peptide_idx), fragment_mz_(fragment_mz)
  {
    follow_up_peaks_ = multiPeak.getFollowUpPeaks();
    //follow_up_peaks_AA_ = multiPeak.getFollowUpPeaksAa();
  }

  MultiFragment::MultiFragment(const OpenMS::MultiFragment& other) = default;

  MultiFragment& MultiFragment::operator=(const OpenMS::MultiFragment& other)
  {
    if(&other == this)
      return *this;
    peptide_idx_ = other.peptide_idx_;
    fragment_mz_ = other.fragment_mz_;
    follow_up_peaks_ = other.follow_up_peaks_;
    //follow_up_peaks_AA_ = other.follow_up_peaks_AA_;
    return *this;
  }

  void MultiFragment::swap(OpenMS::MultiFragment& other)
  {
    follow_up_peaks_.swap(other.follow_up_peaks_);
    //follow_up_peaks_AA_.swap(other.follow_up_peaks_AA_);
    std::swap(peptide_idx_, other.peptide_idx_);
    std::swap(fragment_mz_, other.fragment_mz_);
  }

  size_t MultiFragment::getPeptideIdx() const
  {
    return peptide_idx_;
  }
  double MultiFragment::getFragmentMz() const
  {
    return fragment_mz_;
  }

  const vector<double>& MultiFragment::getFollowUpPeaks() const
  {
    return follow_up_peaks_;
  }
}