// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  Copyright The OpenMS team, Eberhard Karls University Tübingen,
//  ETH Zürich and FU Berlin 2001-2012.
//  This software is released under a BSD license. For a full list of
//  authors, refer to the file AUTHORS. For full licensing conditions
//  refer to the file LICENSE.
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMDECOY_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMDECOY_H

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

namespace OpenMS
{
/**
 @brief This class generates a TargetedExperiment object with decoys based on a TargetedExperiment object

 */

class OPENMS_DLLAPI MRMDecoy: public ProgressLogger
{
public:
  typedef std::vector<OpenMS::TargetedExperiment::Protein> ProteinVectorType;
  typedef std::vector<OpenMS::TargetedExperiment::Peptide> PeptideVectorType;
  typedef std::vector<OpenMS::ReactionMonitoringTransition> TransitionVectorType;

  typedef std::map<String, std::map<String, double> > IonSeries;
  typedef std::map<String, IonSeries> IonSeriesMapType;

  typedef std::map<String, std::vector<const ReactionMonitoringTransition*> > PeptideTransitionMapType;

  MRMDecoy();

  std::pair<String, DoubleReal> getDecoyIon(String ionid,
      std::map<String, std::map<String, DoubleReal> > & decoy_ionseries);

  std::pair<String, double> getTargetIon(double ProductMZ, double mz_threshold,
      std::map<String, std::map<String, double> > target_ionseries);

  inline std::map< String , std::map<String, double> > getIonSeries(
      AASequence sequence, int precursor_charge);

  /// all the following methods are static TODO move them to a procedural Algo part.

  ///TODO add comment
  static std::vector<std::pair<std::string::size_type, std::string> > find_all_tryptic(
      std::string sequence);
  float AASequenceIdentity(const String & sequence, const String & decoy);

  ///AFAIK only the sequence (std::string) of an Peptide object is used here,
  ///TODO pass string only to function
  OpenMS::TargetedExperiment::Peptide shufflePeptide(
      OpenMS::TargetedExperiment::Peptide peptide, double identity_threshold, int seed = -1,
      int maxattempts = 10);

  OpenMS::TargetedExperiment::Peptide reversePeptide(
      OpenMS::TargetedExperiment::Peptide peptide);

  OpenMS::TargetedExperiment::Peptide trypticreversePeptide(
      OpenMS::TargetedExperiment::Peptide peptide);

  int getPeptideItbyId(OpenMS::TargetedExperiment& exp, String id);

  void annotateTransitions(TargetedExperiment &exp, double mz_threshold,
      IonSeriesMapType &IonSeriesMap);

  ///TODO add description
  void restrictTransitions(OpenMS::TargetedExperiment &exp, int min_transitions,
      int max_transitions);

  OpenMS::AASequence getAASequence(const OpenMS::TargetedExperiment::Peptide & peptide);

  //////
  // Non Static functions in this class This seems to be the main entry point //TODO better separate from other code.
  /////

  IonSeriesMapType getIonSeriesMap(TargetedExperiment *exp);

  void generateDecoys(OpenMS::TargetedExperiment& exp,
      OpenMS::TargetedExperiment& dec, String method, String decoy_tag,
      double identity_threshold, double mz_threshold, bool theoretical);
};
}

#endif
