// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/config.h>
#include <iomanip>
#include <iostream>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>

namespace OpenMS
{
  class OPENMS_DLLAPI FLASHTnTFile
{
  public:
    /// write header line for tag output file
    static void writeTagHeader(std::fstream& fs);

    /// write header line for PrSM output file. With @p chimeric (FLASHTnT -chimeric_rounds > 0), the chimeric columns follow the others:
    /// ChimericRound (0: the search of the whole spectrum, r: residual round r), MatchedMasses (deconvolved masses on the PrSM's path),
    /// RemovedMasses (masses removed from the spectrum before the PrSM's round), ResidualNullScore (the best extension score of the other
    /// decoy PrSMs of the same residual search; -1: none or round 0) and ResidualScoreRatio (extension score / ResidualNullScore; -1:
    /// undefined)
    /// With @p targets (FLASHTnT -target), TargetProtein and DecoyMargin follow.
    static void writePrSMHeader(std::fstream& fs, bool chimeric = false, bool targets = false);

    /// write header line for Proteoform output file (@p chimeric: as writePrSMHeader, the values of the representative PrSM)
    static void writeProHeader(std::fstream& fs, bool chimeric = false, bool targets = false);

    /// header for the protein-level output, written before proteoform characterization
    static void writeProteinHeader(std::fstream& fs);

    /// the candidate proteins and their protein-level scores
    static void writeProteins(const std::vector<ProteinHit>& hits, std::fstream& fs);

    /// write the features in regular file output
    static void writeTags(const FLASHTnTAlgorithm& tnt, double flanking_mass_tol, std::fstream& fs);

    static void writePrSMs(const std::vector<ProteinHit>& hits, std::fstream& fs, bool chimeric = false, bool targets = false);

    /**
      @brief the ProForma 2.0 string of a PrSM or proteoform hit; empty if a terminus is undetermined (ProForma cannot express it)

      A variable (searched) modification is written by its name (two names joined by '+' as two tags), a blind shift by its mass: the
      Unimod names listed for a blind shift (ModID) are only those of matching mass. Overlapping ranges, which ProForma does not
      allow, are written as one range with their summed mass.
    */
    static std::string getProFormaString(const ProteinHit& hit);

    /// the fixed modifications of the whole search (FLASHTnT -fixed_mod) as ProForma global modifications, e.g. "<[+57.021464|UNIMOD:4]@C>",
    /// written before every ProForma string (empty: none)
    static void setGlobalFixedModifications(const std::string& proforma_prefix);

    static void writeProteoforms(const std::vector<ProteinHit>& hits, std::fstream& fs, double pro_fdr, bool chimeric = false, bool targets = false);
  };
} // namespace OpenMS