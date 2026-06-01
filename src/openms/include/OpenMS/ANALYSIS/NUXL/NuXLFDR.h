// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLReport.h>
#include <vector>

namespace OpenMS
{

/**
  @brief Class-separated FDR calculation for nucleic-acid cross-links (NuXL).

  Splits PSMs into a "plain peptide" class (@c NuXL:isXL @c == @c 0) and a "cross-link"
  class (@c NuXL:isXL @c != @c 0), computes FDR separately for each, and (optionally)
  writes one filtered idXML per requested XL q-value threshold plus a TSV protein report
  for stringent (@c <= @c 10%) XL FDRs.

  Built on top of @ref FalseDiscoveryRate. Decoys are kept in the FDR result (so downstream
  Percolator and similar tools can still see them); the convenience filter method calls
  @ref IDFilter::removeDecoyHits before writing.

  @ingroup Analysis_ID
*/
class OPENMS_DLLAPI NuXLFDR
{
  public:
    /**
      @brief Construct with a top-hit cap that controls how @ref FalseDiscoveryRate is parametrised.

      @p report_top_hits @c >= @c 2 toggles @c use_all_hits @c = @c true on the underlying
      @ref FalseDiscoveryRate so q-values are computed across all hits per PSM, not only
      the top hit. The value is otherwise stored verbatim and consulted by the
      @ref QValueAtPSMLevel and @ref calculatePeptideAndXLQValueAtPSMLevel methods.

      @param[in] report_top_hits Top-hit cap; values @c >= @c 2 enable the underlying FDR engine's @c use_all_hits mode.
    */
    explicit NuXLFDR(size_t report_top_hits);

    /**
      @brief Partition @p peptide_ids into plain-peptide and cross-link PSM lists.

      For each input @c PeptideIdentification, walks its hits and keeps **only the first**
      hit encountered for each class — even when @c report_top_hits @c >= @c 2 was passed
      to the constructor. The classification is based on the integer meta value
      @c "NuXL:isXL" (zero → plain peptide; anything else → cross-link). An input
      identification contributes to @p pep_pi or @p xl_pi (or both, if it had hits of both
      classes), but never to @p pep_pi twice.

      @p pep_pi and @p xl_pi are cleared on entry.

      @param[in]  peptide_ids Input mixed plain-peptide / cross-link PSMs.
      @param[out] pep_pi      Receives the plain-peptide PSMs (one hit per PSM, the first encountered with @c NuXL:isXL @c == @c 0).
      @param[out] xl_pi       Receives the cross-link PSMs (one hit per PSM, the first encountered with @c NuXL:isXL @c != @c 0).
    */
    void splitIntoPeptidesAndXLs(const PeptideIdentificationList& peptide_ids,
      PeptideIdentificationList& pep_pi,
      PeptideIdentificationList& xl_pi) const;

    /**
      @brief Re-merge a previously split (@ref splitIntoPeptidesAndXLs) pair of lists into one PSM list.

      @p pep_pi entries are appended verbatim. For each @p xl_pi entry, the meta value
      @c "spectrum_reference" is used as the merge key:
        - if no plain-peptide entry shares this spectrum reference, the XL entry is
          appended as a new @c PeptideIdentification;
        - if a plain-peptide entry already exists for this spectrum, the XL hits are
          appended to its hit list and that hit list is re-sorted in place.

      @p peptide_ids is cleared on entry.

      @param[in]  pep_pi      Plain-peptide PSMs (typically from @ref splitIntoPeptidesAndXLs).
      @param[in]  xl_pi       Cross-link PSMs (typically from @ref splitIntoPeptidesAndXLs).
      @param[out] peptide_ids Receives the merged PSM list keyed by @c "spectrum_reference".
    */
    void mergePeptidesAndXLs(const PeptideIdentificationList& pep_pi,
      const PeptideIdentificationList& xl_pi,
      PeptideIdentificationList& peptide_ids) const;

    /**
      @brief Compute PSM-level q-values across @p peptide_ids without splitting by class.

      Wraps @ref FalseDiscoveryRate with @c add_decoy_proteins @c = @c true,
      @c add_decoy_peptides @c = @c true (decoys are retained in the result for downstream
      tools such as Percolator), and @c use_all_hits @c = @c true iff @c report_top_hits_
      passed to the constructor is @c >= @c 2. Also computes the peptide-level q-value
      meta value @c Constants::UserParam::PEPTIDE_Q_VALUE on each hit.

      @param[in,out] peptide_ids PSMs to score in place; each hit receives PSM- and peptide-level q-value meta values.
    */
    void QValueAtPSMLevel(PeptideIdentificationList& peptide_ids) const;

    /**
      @brief Compute PSM- and peptide-level q-values separately for plain and cross-linked PSMs.

      Internally calls @ref splitIntoPeptidesAndXLs and then applies
      @ref FalseDiscoveryRate to each output list with the same parameters as
      @ref QValueAtPSMLevel. The q-values are computed independently on each class so XL
      and non-XL FDR control are not confounded.

      @param[in]  peptide_ids Input mixed plain-peptide / cross-link PSMs.
      @param[out] pep_pi      Plain-peptide PSMs with PSM- and peptide-level q-values annotated.
      @param[out] xl_pi       Cross-link PSMs with PSM- and peptide-level q-values annotated.
    */
    void calculatePeptideAndXLQValueAtPSMLevel(const PeptideIdentificationList& peptide_ids,
      PeptideIdentificationList& pep_pi,
      PeptideIdentificationList& xl_pi) const;

    /**
      @brief Compute separated FDRs, apply class-specific filters, and write per-threshold idXML / TSV reports.

      Pipeline (every step is observable behaviour, not implementation detail):
        -# Calls @ref calculatePeptideAndXLQValueAtPSMLevel to populate @p pep and @p xl_pi
           with FDR-annotated PSMs.
        -# When @p decoy_factor @c != @c 1, divides every cross-link hit's score by
           @p decoy_factor (peptide hits are not touched).
        -# Adds a deterministic score-range-normalised tie-breaker of magnitude
           @c <= @c 1e-5 to every cross-link score so PSMs that originally tied on the
           q-value resolve to a unique ordering. Uses @c "svm_score" when at least one XL
           hit carries it, otherwise @c "NuXL:score". @warning the tie-break divides by
           the @c (max @c - @c min) score range and is undefined when all XL hits share
           the same score.
        -# Calls @ref IDFilter::removeDecoyHits on both @p pep and @p xl_pi.
        -# Filters @p pep on @c PEPTIDE_Q_VALUE when @p peptide_peptide_qvalue_threshold
           is in @c (0, 1); values @c <= @c 0 or @c >= @c 1 disable the filter.
        -# Filters @p pep PSMs by score when @p peptide_PSM_qvalue_threshold is in
           @c (0, 1).
        -# Writes the plain-peptide idXML to @c {out_idxml}{threshold:.4f}_peptides.idXML
           (with @p out_idxml the literal prefix and @c threshold the configured cut
           formatted to four decimal places), with unreferenced proteins pruned.
        -# Normalises @p xl_PSM_qvalue_thresholds: @c 0.0 entries are rewritten to @c 1.0
           (disabled filter = 100% FDR) and the vector is then sorted in @em descending
           order so increasingly stringent filters are applied progressively to the same
           @p xl_pi (filters compound across iterations).
        -# For each (XL-PSM threshold, XL-peptide-level threshold) pair (paired by index;
           sizes must match — otherwise @c Exception::InvalidValue is thrown):
            - filter cross-links on @c PEPTIDE_Q_VALUE when the peptide-level threshold
              is in @c (0, 1);
            - filter cross-links by PSM-level q-value when the XL-PSM threshold is in
              @c (0, 1);
            - prune unreferenced proteins from a fresh copy of @p protein_ids;
            - compute coverage by cross-linked peptides on the first @c ProteinIdentification
              of that copy;
            - write the XL idXML to @c {out_idxml}{threshold:.4f}_XLs.idXML;
            - if the XL FDR is @c <= @c 0.1, also write a TSV protein report to
              @c {out_idxml}proteins{threshold:.4f}_XLs.tsv (skipped for permissive FDRs
              to bound output size).

      @param[in]  protein_ids                       Protein identifications used as the source for the per-file protein lists; unreferenced proteins are pruned per output.
      @param[in]  peptide_ids                       Input PSMs (mixed plain + XL).
      @param[out] pep                               Receives the filtered plain-peptide PSMs (after step 6).
      @param[in]  peptide_PSM_qvalue_threshold      Plain-peptide PSM-level q-value cut; values outside @c (0, 1) disable the filter.
      @param[in]  peptide_peptide_qvalue_threshold  Plain-peptide peptide-level q-value cut; values outside @c (0, 1) disable the filter.
      @param[out] xl_pi                             Receives the filtered cross-link PSMs (state at the most stringent XL threshold).
      @param[in]  xl_PSM_qvalue_thresholds          Cross-link PSM-level q-value thresholds; @c 0.0 is rewritten to @c 1.0 (disabled = 100% FDR). The vector is sorted descending in place.
      @param[in]  xl_peptidelevel_qvalue_thresholds Cross-link peptide-level q-value thresholds paired by index with @p xl_PSM_qvalue_thresholds (must have the same size).
      @param[in]  out_idxml                         Output filename @em prefix; the FDR threshold and a suffix (@c _peptides.idXML / @c _XLs.idXML / @c proteins...XLs.tsv) are appended.
      @param[in]  decoy_factor                      Score scale-down factor for the XL class; ignored when equal to @c 1.

      @throws OpenMS::Exception::InvalidValue when @p xl_PSM_qvalue_thresholds and @p xl_peptidelevel_qvalue_thresholds differ in size.
    */
    void calculatePeptideAndXLQValueAndFilterAtPSMLevel(
      const std::vector<ProteinIdentification>& protein_ids,
      const PeptideIdentificationList& peptide_ids,
      PeptideIdentificationList& pep,
      double peptide_PSM_qvalue_threshold,
      double peptide_peptide_qvalue_threshold,
      PeptideIdentificationList& xl_pi,
      std::vector<double> xl_PSM_qvalue_thresholds,
      std::vector<double> xl_peptidelevel_qvalue_thresholds,
      const String& out_idxml,
      int decoy_factor) const;

  private:
    size_t report_top_hits_; ///< Top-hit cap consulted by the FDR-applying methods to toggle @c use_all_hits.
};

}


