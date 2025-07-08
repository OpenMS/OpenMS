// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLMarkerIonExtractor.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <vector>

namespace OpenMS
{

/// @brief struct to hold a single report line
struct OPENMS_DLLAPI NuXLReportRow
{
  bool no_id;

  // columns
  double rt;
  double original_mz;
  String accessions;
  String peptide;
  String NA;
  Int charge;
  double score;
  int rank;
  double best_localization_score;
  String localization_scores;
  String best_localization;
  double peptide_weight;
  double NA_weight;
  double xl_weight;
  StringList meta_values; // the actual values of exported metadata
  NuXLMarkerIonExtractor::MarkerIonsType marker_ions;
  double abs_prec_error;
  double rel_prec_error;
  double m_H;
  double m_2H;
  double m_3H;
  double m_4H;
  String fragment_annotation;
  String getString(const String& separator) const;
};

/// create header line
struct OPENMS_DLLAPI NuXLReportRowHeader
{
  static String getString(const String& separator, const StringList& meta_values_to_export);
};

/// create PSM report
struct OPENMS_DLLAPI NuXLReport
{
  static std::vector<NuXLReportRow> annotate(
    const PeakMap& spectra, 
    std::vector<PeptideIdentification>& peptide_ids, 
    const StringList& meta_values_to_export,
    double marker_ions_tolerance);
};


/// protein report
struct OPENMS_DLLAPI NuXLProteinReport
{
  static void annotateProteinModificationForTopHits(std::vector<ProteinIdentification>& prot_ids, 
    const std::vector<PeptideIdentification>& peps, 
    TextFile& tsv_file);

  // crosslink efficiency = frequency of the crosslinked amino acid / frequency of the amino acid in all crosslink spectrum matches
  static std::map<char, double> getCrossLinkEfficiency(const std::vector<PeptideIdentification>& peps);

  // returns map of adduct to counts
  static std::map<String, size_t> countAdducts(const std::vector<PeptideIdentification>& peps);

  static void mapAccessionToTDProteins(ProteinIdentification& prot_id, std::map<String, ProteinHit*>& acc2protein_targets, std::map<String, ProteinHit*>& acc2protein_decoys);
};

}
