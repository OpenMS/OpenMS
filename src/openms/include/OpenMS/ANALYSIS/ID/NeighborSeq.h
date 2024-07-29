// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Philipp Wang $
// $Authors: Philipp Wang $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>



namespace OpenMS
{
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
       * @brief Generates a theoretical spectrum for a given peptide sequence.
       * @param peptide_sequence The peptide sequence for which to generate the spectrum.
       * @return The generated theoretical spectrum.
       */
      static MSSpectrum generateSpectrum(const String& peptide_sequence);

      /**
       * @brief Compares two spectra to determine if they share a sufficient number of ions.
       * @param spec1 The first spectrum.
       * @param spec2 The second spectrum.
       * @param min_shared_ion_fraction The tolerance for the proportion of shared ions.
       * @param mz_bin_size The mz_bin_size setting for the comparison.
       * @return True if the spectra share a sufficient number of ions, false otherwise.
       */
      static bool compareSpectra(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& min_shared_ion_fraction, const double& mz_bin_size);

      /**
       * @brief Finds candidate positions based on a given mono-isotopic weight and mass tolerance.
       * @param mass_position_map A map where the key is the mass and the value is a vector of positions.
       * @param mono_weight The mono-isotopic weight to find candidates for.
       * @param mass_tolerance The allowed tolerance for matching the mass.
       * @return A vector of candidate positions within the specified mass tolerance range.
       */
      static std::vector<int> findCandidatePositions(const std::map<double, std::vector<int>>& mass_position_map, const double& mono_weight, const double& mass_tolerance);

      
      /**
       * @brief Creates a map of masses to positions from a vector of peptides. (not protein)
       * @param candidates A vector of AASequence objects.
       * @return A map where the key is the mass and the value is a vector of positions.
       */
      static std::map<double, std::vector<int>> NeighborSeq::createMassPositionMap(const std::vector<AASequence>& candidates);

      /**
       * @brief Finds neighbor peptides that are similar to a given peptide.
       * @param peptide The peptide sequence (from a relevant protein) to find neighbors for.
       * @param neighbor_candidate The list of neighbor sequences to search in.
       * @param mass_tolerance The mass tolerance for neighbor peptides.
       * @param min_shared_ion_fraction The ion tolerance for neighbor peptides.
       * @param mz_bin_size The mz_bin_size setting for the comparison is default 0.05 (the original study suggests ‘high’ (0.05 Da) and ‘low’ (1.0005079
       * Da) mz_bin_size).
       * @return A vector of FASTA entries that are considered neighbor peptides.
       */
      static std::vector<int> findNeighborPeptides(const AASequence& peptides,
                                                   const std::vector<AASequence>& neighbor_candidate,
                                                   const std::vector<int>& candidate_position,
                                                   const double& min_shared_ion_fraction,
                                                   const double& mz_bin_size);

       /**
       * @brief Compares two spectra to determine if they share a sufficient number of ions.
       * @param spec1 The first spectrum.
       * @param spec2 The second spectrum.
       * @param mz_bin_size The mz_bin_size setting for the comparison.
       * @return True if the spectra share a sufficient number of ions, false otherwise.
       */
      static int compareShareSpectra(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& mz_bin_size);
    }; 
}
