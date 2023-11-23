
#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGeneratorNode.h>

#include <OpenMS/DATASTRUCTURES/MultiPeak.h>

#include <utility>
#include <vector>
#include <functional>


namespace OpenMS
{


  class MultiFragment
  {
  public:
    MultiFragment();

    MultiFragment(size_t peptide_idx,
                  double fragment_mz,
                  const std::vector<double>& follow_up);

    MultiFragment(Size peptide_idx, double fragment_mz, const MultiPeak& multiPeak);

    MultiFragment(const MultiFragment& other);

    /// Assignment operator
    MultiFragment& operator=(const MultiFragment& other);

    ///Destructor
    virtual ~MultiFragment() = default;

    /// ValueSwappable
    void swap(MultiFragment& other);


    size_t getPeptideIdx() const;
    double getFragmentMz() const;
    //const std::string& getFollowUpPeaksAa() const;
    const std::vector<double>& getFollowUpPeaks() const;


  protected:
    size_t peptide_idx_;
    double fragment_mz_;
    std::vector<double> follow_up_peaks_;


  };



}
