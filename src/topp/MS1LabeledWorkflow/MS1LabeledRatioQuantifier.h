// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <map>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Channel ratios of MS1-labeled (SILAC, Dimethyl, ...) data, aggregated as medians of ratios

    A labeled experiment measures its channels in one run, so the quantity of interest is their ratio,
    and the ratio of a peptide multiplet is measured far better than either of its channel intensities:
    both channels see the same elution, the same ionization and the same instrument state. This class
    therefore aggregates <em>ratios</em>, the way MaxQuant does, rather than dividing aggregated
    intensities:

      - <b>evidence ratio</b>: per multiplet and run, the intensity of a channel divided by the
        intensity of the reference channel (channel 1, the light one, by default). Only positive
        intensities take part: a channel completed with a zero-intensity dummy feature (the peptide is
        absent) or a not-quantifiable one contributes no ratio, since neither 0 nor infinity is a
        measurement of a ratio.
      - <b>peptide ratio</b>: the median of the evidence ratios of one peptide within one fraction
        group (a fraction group is one labeled sample, measured in one or several fractions).
      - <b>protein group ratio</b>: the median of the peptide ratios of the group's peptides, again
        per fraction group, reported only when at least @p min_ratio_count peptides contribute
        (MaxQuant's "min. ratio count"). The number of contributing peptides is reported alongside.
      - <b>normalized ratio</b>: every ratio divided by the median peptide ratio of its
        (fraction group, channel), i.e. the assumption that most peptides do not change. Medians are
        equivariant under that division, so normalizing the peptide ratios and re-aggregating gives the
        same protein ratios as dividing the protein ratios directly.

    A ratio is therefore <em>not</em> the ratio of the abundances that PeptideAndProteinQuant reports:
    that one is a ratio of per-channel aggregates, which is a different statistic (it weights peptides
    by their intensity, and a single intense peptide can dominate it). Both are useful and both are
    reported by MaxQuant; this class provides the ratio-of-ratios one.

    Ratios are annotated where they belong:
      - on every consensus feature, as the parallel meta values @c evidence_ratio_run (StringList),
        @c evidence_ratio_channel (IntList) and @c evidence_ratio (DoubleList) for its own evidence
        ratios, and @c peptide_ratio_fraction_group, @c peptide_ratio_channel, @c peptide_ratio,
        @c peptide_ratio_normalized and @c peptide_ratio_count for the ratios of its peptide
      - on every indistinguishable protein group, as the parallel data arrays
        @c fraction_group_level_ratio_fraction_group, @c fraction_group_level_ratio_label and
        @c fraction_group_level_ratio_count (integer) plus @c fraction_group_level_ratio and
        @c fraction_group_level_ratio_normalized (float), next to the abundance arrays that
        PeptideAndProteinQuant writes

    Part of the MS1LabeledWorkflow tool rather than of the library: the aggregation rules it implements
    are those of that workflow, and its parameters are the tool's @c ratios section.
  */
  class MS1LabeledRatioQuantifier : public DefaultParamHandler
  {
  public:
    /// One channel ratio of one peptide or protein group, at (fraction group, channel) grain
    struct ChannelRatio
    {
      /// Fraction group the ratio was measured in (1-based)
      unsigned fraction_group = 0;
      /// Channel in the numerator (1-based label of the experimental design; the denominator is the reference channel)
      unsigned channel = 0;
      /// Median of the contributing ratios
      double ratio = 0.0;
      /// @p ratio divided by the median peptide ratio of this (fraction group, channel)
      double normalized_ratio = 0.0;
      /// Number of contributing ratios (evidences for a peptide, peptides for a protein group)
      Size count = 0;
    };

    MS1LabeledRatioQuantifier();

    /**
      @brief Compute the ratios of @p consensus and annotate features and protein groups with them.

      @param[in,out] consensus Linked multiplets, one column per (run, channel); features are annotated with their evidence ratios
      @param[in] design The experimental design; supplies the fraction group of every run
      @param[in,out] proteins The inferred identification run whose indistinguishable groups are annotated
    */
    void run(ConsensusMap& consensus, const ExperimentalDesign& design, ProteinIdentification& proteins);

    /// Peptide ratios of the last run(), by peptide identity
    const std::map<AASequence, std::vector<ChannelRatio>>& getPeptideRatios() const;

    /// Protein group ratios of the last run(), by the group's leading accession
    const std::map<std::string, std::vector<ChannelRatio>>& getProteinGroupRatios() const;

  protected:
    /// Median of @p values (which is sorted in the process); NaN if empty
    static double median_(std::vector<double>& values);

    std::map<AASequence, std::vector<ChannelRatio>> peptide_ratios_;
    std::map<std::string, std::vector<ChannelRatio>> protein_group_ratios_;
  };
} // namespace OpenMS
