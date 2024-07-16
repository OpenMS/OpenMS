// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sven Nahnsen $
// $Authors: Sven Nahnsen, Andreas Bertsch, Chris Bielow, Philipp Wang $
// --------------------------------------------------------------------------

// #define NEIGHBORSEQ_H
#pragma once
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <cmath>
#include <regex>
#include <set>

using namespace OpenMS;

/**
 * @brief The Neighbor Peptide functionality in the TOPPDecoyDatabase class is designed to find peptides from a given set of sequences (FASTA file)
 * that are similar to a target peptide based on mass and spectral characteristics. This section will detail the methods and functionalities
 * specifically related to the Neighbor Peptide search.
 *
 * The paper on subset neighbor search is www.ncbi.nlm.nih.gov/pmc/articles/PMC8489664/
 * PMID: 34236864 PMCID: PMC8489664 DOI: 10.1021/acs.jproteome.1c00483
 */
class OPENMS_DLLAPI NeighborSeq
{

public:
  /**
   * @brief Calculates the monoisotopic mass of a peptide.
   * @param peptide The peptide sequence.
   * @return The monoisotopic mass of the peptide.
   */
  static double calculateMass(const AASequence& peptide);

  /**
   * @brief Generates a theoretical spectrum for a given peptide sequence.
   * @param peptide_sequence The peptide sequence for which to generate the spectrum.
   * @return The generated theoretical spectrum.
   */
  static MSSpectrum generateSpectrum(const String& peptide_sequence);

  /**
   * @brief Compares two spectra to determine if they share a sufficient number of ions.
   * @param spec1 The first spectrum.
   * @param spec2 The second spectrum.
   * @param ion_tolerance The tolerance for the proportion of shared ions.
   * @param resolution The resolution setting for the comparison ("high" or "low").
   * @return True if the spectra share a sufficient number of ions, false otherwise.
   */
  static bool compareSpectra(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& ion_tolerance, const std::string& resolution);

  /**
   * @brief Finds neighbor peptides that are similar to a given peptide.
   * @param peptides The peptide sequence to find neighbors for.
   * @param neighbor_file The list of neighbor sequences to search in.
   * @param mass_tolerance The mass tolerance for neighbor peptides.
   * @param ion_tolerance The ion tolerance for neighbor peptides.
   * @param resolution The resolution setting for the comparison ("high" or "low").
   * @return A vector of FASTA entries that are considered neighbor peptides.
   */
  static std::vector<FASTAFile::FASTAEntry> findNeighborPeptides(const AASequence& peptides,
                                                                 const std::vector<FASTAFile::FASTAEntry>& neighbor_file,
                                                                 const double& mass_tolerance,
                                                                 const double& ion_tolerance,
                                                                 const String& resolution);
};