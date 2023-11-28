
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

    MultiFragment(UInt32 peptide_idx,
                  float fragment_mz,
                  const std::vector<float>& follow_up);

    MultiFragment(UInt32 peptide_idx, float fragment_mz, const MultiPeak& multiPeak);

    MultiFragment(const MultiFragment& other);

    /// Assignment operator
    MultiFragment& operator=(const MultiFragment& other);

    ///Destructor
    virtual ~MultiFragment() = default;

    /// ValueSwappable
    void swap(MultiFragment& other);


    UInt32 getPeptideIdx() const;
    float getFragmentMz() const;
    //const std::string& getFollowUpPeaksAa() const;
    const std::vector<float>& getFollowUpPeaks() const;


  protected:
    UInt32 peptide_idx_;
    float fragment_mz_;
    std::vector<float> follow_up_peaks_;


  };



}
