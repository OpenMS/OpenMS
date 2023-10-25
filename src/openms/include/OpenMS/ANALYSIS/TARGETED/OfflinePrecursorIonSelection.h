// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#pragma once


#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{
  class PeptideIdentification;
  class ProteinIdentification;
  class String;


  /**
    @brief Implements different algorithms for precursor ion selection

    Implements different algorithms for precursor ion selection,
    either based on a whole FeatureMap (e.g. like with LC-MALDI MS data)
    or based on single scans (e.g. with LC-ESI MS data).

    @htmlinclude OpenMS_OfflinePrecursorIonSelection.parameters
  */
  class OPENMS_DLLAPI OfflinePrecursorIonSelection :
    public DefaultParamHandler
  {
public:
    typedef PSLPFormulation::IndexTriple IndexTriple;

    OfflinePrecursorIonSelection();
    ~OfflinePrecursorIonSelection() override;

    /**
      @brief Makes the precursor selection for a given feature map, either feature or scan based.

      @param features Input feature map
      @param experiment Input raw data
      @param ms2 Precursors are added as empty MS2 spectra to this MSExperiment
      @param charges_set Allowed charge states
      @param feature_based If true the selection is feature based, if false it is scan based and the highest signals in each spectrum are chosen
    */
    void makePrecursorSelectionForKnownLCMSMap(const FeatureMap& features,
                                               const PeakMap& experiment,
                                               PeakMap& ms2,
                                               std::set<Int>& charges_set,
                                               bool feature_based);

    /**
      @brief Calculates the mass ranges for each feature and stores them as indices of the raw data.

      @param features Input feature map
      @param experiment Input raw data
      @param indices The boundaries of the features as indices in the raw data
    */
    void getMassRanges(const FeatureMap& features,
                       const PeakMap& experiment,
                       std::vector<std::vector<std::pair<Size, Size> > >& indices);

    void createProteinSequenceBasedLPInclusionList(String include, String rt_model_file, String pt_model_file, FeatureMap& precursors);

    void setLPSolver(LPWrapper::SOLVER solver)
    {
      solver_ = solver;
      std::cout << " LPSolver set to " << solver_ << std::endl;
    }

    LPWrapper::SOLVER getLPSolver()
    {
      return solver_;
    }

private:

    template <typename PairType>
    struct PairComparatorSecondElement
    {
      bool operator()(const PairType & left, const PairType & right) const
      {
        return left.second < right.second;
      }

    };

    typedef std::map<std::pair<double, double>, int, PairComparatorSecondElement<std::pair<double, double> > > ExclusionListType_;

    /**
      @brief Calculate the sum of intensities of relevant features for each scan separately.
    */
    void calculateXICs_(const FeatureMap& features,
                        const std::vector<std::vector<std::pair<Size, Size> > >& mass_ranges,
                        const PeakMap& experiment,
                        const std::set<Int>& charges_set,
                        std::vector<std::vector<std::pair<Size, double> > >& xics);

    /**
      @brief Eliminates overlapping peaks.
    */
    void checkMassRanges_(std::vector<std::vector<std::pair<Size, Size> > >& mass_ranges,
                          const PeakMap& experiment);

    /// reduce scan count for each entry, and remove every entry which has reached 0 counts
    void updateExclusionList_(ExclusionListType_& exclusion_list) const;

    LPWrapper::SOLVER solver_;
  };
 
}


