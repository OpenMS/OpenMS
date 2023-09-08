// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{

/// @brief struct to hold a single report line
struct OPENMS_DLLAPI RNPxlReportRow
{
  bool no_id;
  double rt;
  double original_mz;
  String accessions;
  String RNA;
  String peptide;
  double best_localization_score;
  String localization_scores;
  String best_localization;
  Int charge;
  double score;
  double peptide_weight;
  double RNA_weight;
  double xl_weight;
  double abs_prec_error;
  double rel_prec_error;
  RNPxlMarkerIonExtractor::MarkerIonsType marker_ions;
  double m_H;
  double m_2H;
  double m_3H;
  double m_4H;  
  int rank;
  String getString(const String& separator) const;

};

/// create header line
struct OPENMS_DLLAPI RNPxlReportRowHeader
{
  static String getString(const String& separator);
};

/// create report
struct OPENMS_DLLAPI RNPxlReport
{
  static std::vector<RNPxlReportRow> annotate(const PeakMap& spectra, std::vector<PeptideIdentification>& peptide_ids, double marker_ions_tolerance);
};

}


