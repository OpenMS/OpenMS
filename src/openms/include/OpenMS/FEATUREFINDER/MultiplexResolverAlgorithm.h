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
#include <unordered_map>
#include <vector>

namespace OpenMS
{
  class AASequence;
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

    With @c algorithm:use_blacklist_owners, and a blacklist that records who reserved each peak (hasBlacklistOwners(),
    as FeatureFinderMultiplexAlgorithm writes it), only peaks reserved by another multiplet make a channel not
    quantifiable (isForeign()): peaks of the multiplet's own channels (their isotopes) and of no reported feature (a
    pattern the detection dropped later) do not. A dummy feature of intensity 0 then means that no other multiplet
    occupies its position, not that the peptide is absent, so the caller has to measure it (MS1LabeledWorkflow
    re-quantifies such channels). All blacklisted peaks within the m/z tolerance count, on either side of the
    expected positions, and channels that MS1LabeledWorkflow added for an identification (intensity 0, no unique id)
    are checked like dummy features. Otherwise (the default) the owners are ignored: a dummy feature is not
    quantifiable if the first blacklisted peak at or above one of its expected positions is within the tolerance,
    and added channels are not checked.

    Only the first peptide identification of a consensus feature is taken into account, so the map should
    have been reduced to one identification per feature (IDConflictResolverAlgorithm::resolve()) first.
    checkIdentifications() tells beforehand which of several identifications agree with a multiplet, so that
    the reduction can prefer one of those over a better-scoring one that would make the multiplet a conflict.
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
                 If its owners are used (@c algorithm:use_blacklist_owners and hasBlacklistOwners()), the feature handles of
                 @p map_in must carry the unique ids of their features, as FeatureFinderMultiplexAlgorithm::getConsensusMap() does.

      @throws Exception::MissingInformation if an identified feature lacks the @c map_index meta value
      @throws Exception::InvalidValue if the @c map_index does not name a feature handle of the consensus feature
    */
    void resolve(const ConsensusMap& map_in, ConsensusMap& map_out, ConsensusMap& map_conflicts,
                 const MSExperiment& blacklist = MSExperiment()) const;

    /**
      @brief Which peptide identifications of @p consensus agree with its detected mass shifts?

      This is the check resolve() applies to the (first) identification of a multiplet, for every
      identification: the label set of its first peptide hit has to occur in a theoretical pattern whose
      mass shifts explain every detected feature, relative to the feature named by the identification's
      @c map_index. Use it to choose, among several identifications mapped onto one multiplet, one that
      does not turn the multiplet into a conflict.

      @return One flag per peptide identification of @p consensus, in order. Identifications without
              hits, without @c map_index, or whose @c map_index names no feature handle are not consistent.
    */
    std::vector<bool> checkIdentifications(const ConsensusFeature& consensus) const;

    /**
      @brief Where the peptide of @p sequence appears in every channel.

      The theoretical pattern of the peptide is the one whose label set at some position equals the labels
      of @p sequence and whose channels label exactly the residues it has (see matchesResidues_()). Its
      mass shifts, relative to the position of @p sequence, give the m/z of the peptide at charge
      @p charge in every channel.

      @param[in] sequence Identified peptidoform, with its labels
      @param[in] charge Charge state (non-zero)
      @param[out] position Position (channel) of @p sequence itself within the returned m/z values
      @return The m/z per channel (pattern position), empty if no theoretical pattern fits @p sequence
    */
    std::vector<double> expectedChannelMZs(const AASequence& sequence, int charge, Size& position) const;

    /**
      @name Owners of blacklisted peaks

      FeatureFinderMultiplexAlgorithm records which feature reserved each blacklisted peak, as two integer data
      arrays aligned with the peaks of each blacklist spectrum: @c owner_unique_id_high and @c owner_unique_id_low,
      the upper and lower 32 bits of the feature's unique id (an integer data array holds 32-bit values). 0 in
      both means that no reported feature owns the peak (its multiplet was dropped after the peak was reserved).
      The arrays survive mzML files (as non-standard data arrays).
    */
    ///@{
    /// Does every spectrum of @p blacklist that has peaks record their owners (and is there a peak)?
    static bool hasBlacklistOwners(const MSExperiment& blacklist);

    /// The unique id of the feature that reserved peak @p peak_index of @p spectrum (UniqueIdInterface::INVALID if none, or if not recorded)
    static UInt64 blacklistOwner(const MSSpectrum& spectrum, Size peak_index);

    /// Record @p owners (one unique id per peak, INVALID if none) as the owners of the peaks of @p spectrum
    static void setBlacklistOwners(MSSpectrum& spectrum, const std::vector<UInt64>& owners);

    /// How the owner of a blacklisted peak relates to a multiplet
    enum class OwnerRelation
    {
      ORPHAN,        ///< no reported feature owns the peak
      OWN,           ///< a detected channel of the multiplet itself
      SAME_PEPTIDE,  ///< a channel of another multiplet identified as the same peptide (and charge)
      UNIDENTIFIED,  ///< a channel of an unidentified multiplet (or of none in the map)
      OTHER_PEPTIDE  ///< a channel of a multiplet identified as another peptide
    };

    /// Peptide (as compared by classifyBlacklistOwner()) of each detected channel of the identified multiplets of a map, by unique id
    using OwnerPeptides = std::unordered_map<UInt64, std::string>;

    /**
      @brief The peptide of every detected channel (feature handle with a unique id) of the identified multiplets in @p multiplets

      A multiplet's peptide is that of its first identification: the sequence without the labels of section
      @c algorithm:labels (other modifications kept) and the charge of the multiplet. So the label states of a
      peptide share it, while another sequence of the same composition or another modification state does not.
    */
    OwnerPeptides ownerPeptides(const ConsensusMap& multiplets) const;

    /// The peptide of @p hit at charge @p charge, as ownerPeptides() names it (the sequence without the labels, and the charge)
    std::string peptideKey(const PeptideHit& hit, int charge) const;

    /**
      @brief How the feature @p owner, which reserved a blacklisted peak, relates to the multiplet @p consensus

      @param[in] owner Unique id of the owner (blacklistOwner())
      @param[in] consensus The multiplet; its detected channels carry the unique ids of their features
      @param[in] peptides ownerPeptides() of the map the multiplets were detected in (incl. @p consensus)
    */
    static OwnerRelation classifyBlacklistOwner(UInt64 owner, const ConsensusFeature& consensus, const OwnerPeptides& peptides);

    /**
      @brief How the feature @p owner, which reserved a blacklisted peak, relates to the peptide @p peptide (of an identification without a multiplet)

      As above, for a peptide that has no multiplet of its own (never OWN).

      @param[in] owner Unique id of the owner (blacklistOwner())
      @param[in] peptide The peptide, as peptideKey() names it
      @param[in] peptides ownerPeptides() of the map the multiplets were detected in
    */
    static OwnerRelation classifyBlacklistOwner(UInt64 owner, const std::string& peptide, const OwnerPeptides& peptides);

    /**
      @brief Does a peak reserved by an owner of this relation belong to another multiplet? (UNIDENTIFIED, SAME_PEPTIDE, OTHER_PEPTIDE)

      Such a peak makes a missing channel at its position not quantifiable, and a channel must not use it. Another multiplet of
      the same peptide (e.g. its missing channel, detected on its own) is quantified itself: measuring the channel from its peaks
      would count them twice.
    */
    static bool isForeign(OwnerRelation relation);
    ///@}

  protected:
    /// The generator of the theoretical mass-shift patterns, as configured in the sections @c algorithm and @c labels
    MultiplexDeltaMassesGenerator createGenerator_() const;

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
      @param[in] detected_delta_mass_at_label_set Detected mass shift of the identified feature (relative to the first feature)
      @param[in] pattern Theoretical pattern
      @param[in] theoretical_delta_mass_at_label_set Theoretical mass shift at which the label set was matched
      @param[out] delta_mass_matched Which delta masses of the theoretical pattern were matched
      @return true if every detected delta mass has a theoretical counterpart within the mass tolerance
    */
    bool matchDeltaMasses_(const ConsensusFeature& consensus, double detected_delta_mass_at_label_set,
                           const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                           double theoretical_delta_mass_at_label_set, std::vector<bool>& delta_mass_matched) const;

    /**
      @brief Find the theoretical delta mass pattern matching the detected one

      @param[in] consensus Detected pattern
      @param[in] detected_delta_mass_at_label_set Detected mass shift of the identified feature (relative to the first feature)
      @param[in] sequence Sequence annotated to the detected pattern
      @param[in] label_set Label set extracted from @p sequence
      @param[in] theoretical_patterns All theoretical patterns
      @param[in] samples Labels of every sample (MultiplexDeltaMassesGenerator::getSamplesLabelsList())
      @param[out] delta_mass_matched Which delta masses of the matching pattern were matched
      @param[out] index_label_set Index within the matching pattern at which the label set was matched
      @return Index of the matching pattern, or -1
    */
    int findMatchingPattern_(const ConsensusFeature& consensus, double detected_delta_mass_at_label_set,
                             const AASequence& sequence, const MultiplexDeltaMasses::LabelSet& label_set,
                             const std::vector<MultiplexDeltaMasses>& theoretical_patterns,
                             const std::vector<std::vector<std::string>>& samples,
                             std::vector<bool>& delta_mass_matched, int& index_label_set) const;

    /**
      @brief Does every channel of @p pattern label exactly the residues @p sequence has?

      Metabolic labels (SILAC Lys/Arg, Leu3) sit on every residue of their kind, so in a pattern that
      fits the peptide, each channel carries one Lys label per K, one Arg label per R and one Leu label
      per L, for every residue kind its sample labels - none if the peptide has no residue of that kind.
      The sample of a channel is the one in @p samples whose labels include all labels of the channel.
      E.g. with samples [][Lys8,Arg10], a peptide with one K and one R fits only the pattern whose heavy
      channel carries Lys8 and Arg10, not the one carrying Lys8 alone. This tells the pattern of an
      unlabelled identification, whose label set ("no_label") occurs in every pattern. Unlabelled channels
      impose nothing, and other labels (Dimethyl, ICPL, plain mass shifts) are not counted.

      @param[in] pattern Theoretical pattern
      @param[in] sequence Identified peptidoform
      @param[in] samples Labels of every sample (MultiplexDeltaMassesGenerator::getSamplesLabelsList())
    */
    static bool matchesResidues_(const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern, const AASequence& sequence,
                                 const std::vector<std::vector<std::string>>& samples);

    /// m/z of the lightest peptide of the complete multiplet, given the incomplete one at @p mz
    static double findNewMZ_(double mz, int charge, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                             const std::vector<bool>& delta_mass_matched);

    /// Is one of the first three isotopes of a (dummy) feature at this position blacklisted? (anonymous blacklist: the first peak at or above each position)
    bool isBlacklisted_(const MSExperiment& blacklist, double rt, double mz, size_t charge, const ConsensusFeature& consensus) const;

    /// What the owner-aware blacklist check found in resolve(), for its log
    struct OwnerCheckCounts_
    {
      Size dummies = 0;        ///< dummy features checked
      Size added = 0;          ///< channels added for an identification checked
      Size dummies_vetoed = 0; ///< dummy features found not quantifiable
      Size added_vetoed = 0;   ///< added channels found not quantifiable
      Size other_peptide = 0;  ///< not quantifiable: a peak reserved by another identified peptide
      Size same_peptide = 0;   ///< not quantifiable: a peak reserved by another multiplet of the same peptide (no other peptide)
      Size unidentified = 0;   ///< not quantifiable: peaks reserved by unidentified multiplets only
      Size usable = 0;         ///< reserved peaks, but none foreign
      Size own = 0;            ///< ... of these, with a peak of the multiplet's own channels
      Size orphan = 0;         ///< ... with a peak of no reported feature
      Size anonymous = 0;      ///< dummy features the check without owners finds blacklisted
    };

    /**
      @brief Is a peak reserved by another multiplet within the m/z tolerance of one of the first three isotopes at this position?

      The owner-aware check of a (dummy or added) channel of @p consensus: every blacklisted peak in the RT and m/z
      tolerance counts (both sides), and only foreign owners (isForeign()) make the channel not quantifiable.
    */
    bool isReservedByOther_(const MSExperiment& blacklist, double rt, double mz, size_t charge, const ConsensusFeature& consensus,
                            const OwnerPeptides& peptides, OwnerCheckCounts_& counts) const;

    /**
      @brief Complete an incomplete multiplet with dummy features

      @param[in] consensus Possibly incomplete multiplet
      @param[in] pattern Matching theoretical pattern
      @param[in] delta_mass_matched Which delta masses of the pattern were matched
      @param[in] index_label_set Index within the pattern at which the label set was matched
      @param[in] blacklist Blacklisted spectral peaks (may be empty)
      @param[in] peptides If not null, the owners of the blacklisted peaks are used and these are the ownerPeptides() of the
                 map: owner-aware check of dummy features and of channels added for an identification
      @param[in,out] counts Counts of the owner-aware check (used only with @p peptides)
      @return The complete multiplet
    */
    ConsensusFeature completeConsensus_(const ConsensusFeature& consensus, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                                        const std::vector<bool>& delta_mass_matched, int index_label_set,
                                        const MSExperiment& blacklist, const OwnerPeptides* peptides = nullptr,
                                        OwnerCheckCounts_* counts = nullptr) const;

    /// Mass tolerance in Da for matching detected to theoretical mass shifts
    double mass_tolerance_ = 0.1;
    /// m/z tolerance in ppm for the blacklist check
    double mz_tolerance_ = 10.0;
    /// RT tolerance in seconds for the blacklist check
    double rt_tolerance_ = 5.0;
    /// Use the owners of blacklisted peaks, if recorded ('algorithm:use_blacklist_owners')
    bool use_blacklist_owners_ = false;

    void updateMembers_() override;
  };
} // namespace OpenMS
