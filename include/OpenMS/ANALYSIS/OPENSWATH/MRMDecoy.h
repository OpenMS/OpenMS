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

#include <map>
#include <string>
#include <vector>
#include <utility>  // for pair

namespace OpenMS
{
/**
 @brief This class generates a TargetedExperiment object with decoys based on a TargetedExperiment object

 This class implements functions to generate decoy peptides with decoy
 transitions from a set of real peptides with real transitions. For this, several 

 */

class OPENMS_DLLAPI MRMDecoy: public ProgressLogger
{
 public:
  MRMDecoy() {} // empty, no members

  /**
    @brief Generate decoys from an TargetedExperiment

    Will generate decoy peptides for each target peptide provided in exp and
    write them into the decoy experiment. 

    Valid methods: reverse, trypticreverse, shuffle

    If theoretical is true, the target transitions will be returned but their
    masses will be adjusted to match the theoretical value of the fragment ion
    that is the most likely explanation for the product.

    mz_threshold is used for the matching of theoretical ion series to the observed one

  */
  void generateDecoys(OpenMS::TargetedExperiment& exp,
      OpenMS::TargetedExperiment& dec, String method, String decoy_tag,
      double identity_threshold, double mz_threshold, bool theoretical);

  /**
    @brief Remove transitions s.t. all peptides have a defined set of transitions.

    All transitions of a peptide above max_transitions get deleted, all
    peptides with less than min_transitions also get deleted.

  */
  void restrictTransitions(OpenMS::TargetedExperiment &exp, int min_transitions,
      int max_transitions);

 private:
  typedef std::vector<OpenMS::TargetedExperiment::Protein> ProteinVectorType;
  typedef std::vector<OpenMS::TargetedExperiment::Peptide> PeptideVectorType;
  typedef std::vector<OpenMS::ReactionMonitoringTransition> TransitionVectorType;

  typedef std::map<String, std::map<String, double> > IonSeries;
  typedef std::map<String, IonSeries> IonSeriesMapType;

  typedef std::map<String, std::vector<const ReactionMonitoringTransition*> > PeptideTransitionMapType;

  /**
    @brief Selects a decoy ion from a set of ions.
  */
  std::pair<String, DoubleReal> getDecoyIon(String ionid,
      std::map<String, std::map<String, DoubleReal> > & decoy_ionseries);

  /**
    @brief Selects a target ion from a set of ions.
  */
  std::pair<String, double> getTargetIon(double ProductMZ, double mz_threshold,
      std::map<String, std::map<String, double> > target_ionseries);

  /**
    @brief Generate all ion series for an input AASequence

    Currently generated are:

    bionseries, bionseries_isotopes, bionseries_loss, bionseries_isotopes_loss, 
    yionseries, yionseries_isotopes, yionseries_loss, yionseries_isotopes_loss, 
    aionseries, aionseries_isotopes
  */
  std::map< String , std::map<String, double> > getIonSeries(
      AASequence sequence, int precursor_charge);

  /**
    @brief Find all tryptic sites in a sequence
  */
  std::vector<std::pair<std::string::size_type, std::string> > find_all_tryptic(
      std::string sequence);

  /**
    @brief Compute (Manhatten?) distance between two sequences
  */
  float AASequenceIdentity(const String & sequence, const String & decoy);

  /**
    @brief Shuffle a peptide (with its modifications) sequence

    This function will shuffle the given peptide sequences and its
    modifications such that the resulting sequence similarity is below
    identity_threshold (TODO: how is seq. similiary measured?). 

    TODO : what is the failure scenario if it cannot be achieved in maxattempts?
  */
  OpenMS::TargetedExperiment::Peptide shufflePeptide(
      OpenMS::TargetedExperiment::Peptide peptide, double identity_threshold, int seed = -1,
      int maxattempts = 10);

  /**
    @brief Pseudo-reverse a peptide sequence (with its modifications)

    Pseudo reverses a peptide sequence, leaving the last AA constant
  */
  OpenMS::TargetedExperiment::Peptide reversePeptide(
      OpenMS::TargetedExperiment::Peptide peptide);

  /**
    @brief Reverse a peptide sequence (with its modifications)
  */
  OpenMS::TargetedExperiment::Peptide trypticreversePeptide(
      OpenMS::TargetedExperiment::Peptide peptide);

  /**
    @brief get AASequence from a peptide
  */
  OpenMS::AASequence getAASequence(const OpenMS::TargetedExperiment::Peptide & peptide);
};
}

#endif
