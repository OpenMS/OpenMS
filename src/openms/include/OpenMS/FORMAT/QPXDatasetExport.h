// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>

#include <string>
#include <vector>

namespace OpenMS
{
  class PeptideIdentificationList;

/**
  @brief The parts of a QPX dataset that a search engine adds to the views of QPXFile

  A tool that writes a QPX dataset (e.g. ProSE with @c -out_qpx) takes the mandatory columns of the psm view from
  QPXFile and adds the optional ones of this class with ArrowIOHelpers::withColumns(). Columns are
  ArrowIOHelpers::ParquetColumn, so the caller needs no Arrow headers.

  All members are static.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI QPXDatasetExport
{
public:
  /// A column of a QPX view, without Arrow types
  using Column = ArrowIOHelpers::ParquetColumn;

  /**
    @brief The name of a score or meta value in the snake_case of the QPX @c additional_scores

    Percolator's scores get their names of bigbio/qpx#355 (e.g. "MS:1001492" -> "percolator_score"), "q-value" is
    "global_qvalue", others are converted ("ln(hyperscore)" -> "ln_hyperscore", "DeltaMass" -> "delta_mass",
    "COMET:xcorr" -> "comet_xcorr").

    @param[in] name The score type or meta value key
    @return the QPX name; "" for the posterior error probability, which has a column of its own
  */
  static std::string scoreName(const std::string& name);

  /**
    @brief The optional columns of the QPX psm view, one entry per row of the view (the best hit of each element of
    @p peptide_ids with hits; elements without hits have no row)

    The columns, in this order:
      - @c additional_scores (SCORE_LIST): the main score and every numeric meta value of the hit, named by
        scoreName() (each name once, the posterior error probability and non-score meta values such as
        @c target_decoy left out), with @c higher_better for the scores whose direction is known;
      - @c spectrum_reference (UTF8): the native ID; null if the element has no spectrum (@p has_spectrum) or no
        native ID;
      - @c isotope_error (INT16): the meta value @c isotope_error (or @c IsotopeError) as observed minus theoretical
        (also so in @c additional_scores); null without one;
      - @c protein_mappings (STRUCT_LIST of @c accession, @c start, @c end, @c aa_before, @c aa_after): one entry per
        peptide evidence, 1-based positions (null if unknown), "-" for a protein terminus, null for an unknown residue.

    @param[in] peptide_ids The identifications of one run, in the order of the rows of the psm view
    @param[in] has_spectrum Per element of @p peptide_ids whether the mz view has its spectrum; empty: all true
    @param[in] isotope_observed_minus_theoretical Whether the isotope error meta values are observed minus
               theoretical (false: theoretical minus observed, their sign is flipped)
    @return the four columns
    @throws Exception::IllegalArgument if @p has_spectrum is neither empty nor of the size of @p peptide_ids
  */
  static std::vector<Column> psmColumns(const PeptideIdentificationList& peptide_ids,
                                        const std::vector<bool>& has_spectrum,
                                        bool isotope_observed_minus_theoretical);
};

} // namespace OpenMS
