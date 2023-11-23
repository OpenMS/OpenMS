#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGeneratorNode.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/MultiPeak.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <utility>
#include <vector>


using namespace std;

namespace OpenMS
{

  MultiPeak::MultiPeak() : peak_(), score_(0), follow_up_peaks({})
  {
  }
  MultiPeak::MultiPeak(OpenMS::Peak1D peak, double score) : peak_(peak), score_(score), follow_up_peaks()
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
  void MultiPeak::addScore(double score)
  {
    score_ += score;
  }

  const Peak1D& MultiPeak::getPeak() const
  {
    return peak_;
  }
  double MultiPeak::getScore() const
  {
    return score_;
  }
  const string& MultiPeak::getFollowUpPeaksAa() const
  {
    return follow_up_peaks_AA;
  }
  const vector<double>& MultiPeak::getFollowUpPeaks() const
  {
    return follow_up_peaks;
  }

  void MultiPeak::addFollowUpPeak(double distance, const std::string &AA)
  {
    follow_up_peaks_AA += AA;
    follow_up_peaks.push_back(distance);
  }
}