// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Lars Nilse, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <map>
#include <string>
#include <vector>

namespace OpenMS
{
  class MultiplexDeltaMassesGenerator;

  /**
    @brief Completes peptide multiplets and resolves quant/ID conflicts within them.

    Tools such as FeatureFinderMultiplexAlgorithm detect peptide feature multiplets (pairs, triplets, ...) in
    MS1-labeled data (SILAC, Dimethyl, ...). Once the multiplets carry peptide sequences (e.g. via IDMapper),
    this class consolidates quantitative and sequence information in two steps:

    - Multiplets whose observed mass shifts contradict the labels found in the annotated sequence are
      moved to a conflict map. Example: a sequence carrying a single Arg10 label mapped to the light
      feature of a SILAC triplet - either the feature detection or the sequence is wrong.
    - Incomplete multiplets (e.g. only the heavy partner of a pair was detected) are completed with
      dummy features. A dummy feature gets intensity 0 when nothing was blacklisted around its
      position during feature detection (the peptide is absent), and NaN (not quantifiable) when
      the region was blacklisted, i.e. another feature overlaps with it.

    Only the first peptide identification of a consensus feature is taken into account, so the map should
    have been reduced to one identification per feature (IDConflictResolverAlgorithm::resolve()) first.
    Multiplets without sequence annotation are written to the conflict map unchanged.

    The identifications must carry the map index of the feature handle they were mapped to in the meta
    value @c map_index (IDMapper with @c annotate_ids_with_subelements). When a multiplet is completed,
    the new map index of the identified feature is recorded in the meta value @c map_index of the first
    peptide hit.

    Parameters: section @c algorithm holds the label specification and tolerances, section @c labels the
    mass shift of every known label (see MultiplexDeltaMassesGenerator).

    @htmlinclude OpenMS_MultiplexResolverAlgorithm.parameters

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI MultiplexResolverAlgorithm : public DefaultParamHandler
  {
  public:
    /// Constructor; registers the @c algorithm and @c labels parameter sections
    MultiplexResolverAlgorithm();

    /**
      @brief Split the annotated multiplets of @p map_in into resolved and conflicting ones.

      @p map_out and @p map_conflicts inherit the meta data of @p map_in (column headers, protein
      identifications, unassigned peptide identifications, data processing) and receive the
      consensus features that could resp. could not be matched to a theoretical mass-shift pattern.
      Both output maps get fresh unique ids for their features, and the column header sizes of
      @p map_out are set to its number of features.

      @param[in] map_in Multiplets with (at most one) peptide identification each
      @param[out] map_out Complete, consistent multiplets; cleared first
      @param[out] map_conflicts Multiplets without identification or with conflicting quant/ID information; cleared first
      @param[in] blacklist Spectral peaks blacklisted during feature detection (FeatureFinderMultiplexAlgorithm::getBlacklist()).
                 Used to decide whether a dummy feature is absent (intensity 0) or not quantifiable (NaN). May be empty.

      @throws Exception::MissingInformation if an identified feature lacks the @c map_index meta value
      @throws Exception::InvalidValue if the @c map_index does not name a feature handle of the consensus feature
    */
    void resolve(const ConsensusMap& map_in, ConsensusMap& map_out, ConsensusMap& map_conflicts,
                 const MSExperiment& blacklist = MSExperiment()) const;

  protected:
    /// Relative delta mass between the first feature handle and the handle with map index @p idx (NaN if absent)
    static double deltaMassFromMapIndex_(const ConsensusFeature::HandleSetType& feature_handles, unsigned idx);

    /**
      @brief Mass shift of the theoretical pattern whose label set equals @p label_set

      @param[in] pattern Theoretical pattern
      @param[in] label_set Label set of the detected pattern
      @param[out] index_label_set Index within the pattern at which the label set was matched
      @return Mass shift at the match, or NaN if no delta mass carries this label set
    */
    static double matchLabelSet_(const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                                 const MultiplexDeltaMasses::LabelSet& label_set, int& index_label_set);

    /**
      @brief Check whether all delta masses of the detected pattern match one of the theoretical pattern

      @param[in] consensus Detected pattern
      @param[in] pattern Theoretical pattern
      @param[in] theoretical_delta_mass_at_label_set Theoretical mass shift at which the label set was matched
      @param[out] delta_mass_matched Which delta masses of the theoretical pattern were matched
      @return true if every detected delta mass has a theoretical counterpart within the mass tolerance
    */
    bool matchDeltaMasses_(const ConsensusFeature& consensus, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                           double theoretical_delta_mass_at_label_set, std::vector<bool>& delta_mass_matched) const;

    /**
      @brief Find the theoretical delta mass pattern matching the detected one

      @param[in] consensus Detected pattern
      @param[in] label_set Label set extracted from the sequence annotated to the detected pattern
      @param[in] theoretical_patterns All theoretical patterns
      @param[out] delta_mass_matched Which delta masses of the matching pattern were matched
      @param[out] index_label_set Index within the matching pattern at which the label set was matched
      @return Index of the matching pattern, or -1
    */
    int findMatchingPattern_(const ConsensusFeature& consensus, const MultiplexDeltaMasses::LabelSet& label_set,
                             const std::vector<MultiplexDeltaMasses>& theoretical_patterns,
                             std::vector<bool>& delta_mass_matched, int& index_label_set) const;

    /// m/z of the lightest peptide of the complete multiplet, given the incomplete one at @p mz
    static double findNewMZ_(double mz, int charge, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                             const std::vector<bool>& delta_mass_matched);

    /// Is one of the first three isotopes of a (dummy) feature at this position blacklisted?
    bool isBlacklisted_(const MSExperiment& blacklist, double rt, double mz, size_t charge) const;

    /**
      @brief Complete an incomplete multiplet with dummy features

      @param[in] consensus Possibly incomplete multiplet
      @param[in] pattern Matching theoretical pattern
      @param[in] delta_mass_matched Which delta masses of the pattern were matched
      @param[in] index_label_set Index within the pattern at which the label set was matched
      @param[in] blacklist Blacklisted spectral peaks (may be empty)
      @return The complete multiplet
    */
    ConsensusFeature completeConsensus_(const ConsensusFeature& consensus, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                                        const std::vector<bool>& delta_mass_matched, int index_label_set,
                                        const MSExperiment& blacklist) const;

    /// Mass tolerance in Da for matching detected to theoretical mass shifts
    double mass_tolerance_ = 0.1;
    /// m/z tolerance in ppm for the blacklist check
    double mz_tolerance_ = 10.0;
    /// RT tolerance in seconds for the blacklist check
    double rt_tolerance_ = 5.0;

    void updateMembers_() override;
  };
} // namespace OpenMS
