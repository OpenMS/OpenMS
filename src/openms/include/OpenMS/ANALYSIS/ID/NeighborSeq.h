// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Philipp Wang $
// $Authors: Chris Bielow, Philipp Wang $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <vector>
#include <map>



namespace OpenMS
{
    /**
       @brief The Neighbor Peptide functionality is designed to find peptides (neighbors) in a given set of sequences (FASTA file) that are
              similar to a target peptide (aka relevant peptide) based on mass and spectral characteristics. This provides more power
              when searching complex samples, when only a subset of the peptides/proteins is of interest.
     
       The paper on subset neighbor search is www.ncbi.nlm.nih.gov/pmc/articles/PMC8489664/
       DOI: 10.1021/acs.jproteome.1c00483
     */
    class OPENMS_DLLAPI NeighborSeq
    {

    public:
      /// Constructor
      /// @param digested_relevant_peptides A vector of digested relevant peptides
      NeighborSeq(std::vector<AASequence>&& digested_relevant_peptides);

      /**
       * @brief Generates a theoretical spectrum for a given peptide sequence with b/y ions at charge 1.
       * 
       * Includes all b and y ions with charge 1 (even the prefix ions, e.g. b1), but no internal ions.
       * 
       * @param peptide_sequence The peptide sequence for which to generate the spectrum.
       * @return The generated theoretical spectrum.
       */
      MSSpectrum generateSpectrum(const AASequence& peptide_sequence);

      /**
       * @brief Compares two spectra to determine if they share a sufficient number of ions.
       *
       * All peaks are considered. Use generateSpectrum() to generate theoretical spectra with b/y ions.
       *
       * @param spec1 The first theoretical spectrum.
       * @param spec2 The second theoretical spectrum.
       * @param min_shared_ion_fraction The minimal required proportion of shared ions in [0, 1]
       * @param mz_bin_size Bin size for the m/z values, which determines if two peaks are considered to be the same (typically, 0.05 for high resolution and 1.0005079 for low resolution).
       * @return True if the spectra share a sufficient number of ions, false otherwise.
       */
      static bool isNeighborSpectrum(const MSSpectrum& spec1, const MSSpectrum& spec2, const double min_shared_ion_fraction, const double mz_bin_size);
      /**
       * @brief Compute the number of shared ions between two spectra 
       *
       * All peaks are considered. Use generateSpectrum() to generate theoretical spectra with b/y ions.
       * 
       * @param spec1 The first theoretical spectrum.
       * @param spec2 The second theoretical spectrum.
       * @param mz_bin_size Bin size for the m/z values, which determines if two peaks are considered to be the same.
       * @return The number of shared ions
       */
      static int computeSharedIonCount(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& mz_bin_size);

      /**
       * @brief Is this peptide a neighbor to one of the relevant peptides?
       * 
       * Also updates the internal statistics, which can be retrieved using getNeighborStats().
       * 
       * @param neighbor_candidate The peptide sequence (from a neighbor protein) to compare against the internal relevant peptides (see constructor).
       * @param mass_tolerance_pc Maximal precursor mass difference (in Da or ppm; see 'mass_tolerance_pc_ppm') between neighbor and relevant peptide.
       * @param mass_tolerance_pc_ppm Is 'mass_tolerance_pc' in Da or ppm?
       * @param min_shared_ion_fraction The ion tolerance for neighbor peptides.
       * @param mz_bin_size Bin size for spectra m/z comparison (the original study suggests 0.05 Th for high-res and 1.0005079 Th for low-res spectra).
       * @return true if @p neighbor_candidate is neighbor to one or more relevant peptides, false otherwise.
       */
      bool isNeighborPeptide(const AASequence& neighbor_candidate,
                             const double mass_tolerance_pc,
                             const bool mass_tolerance_pc_ppm,
                             const double min_shared_ion_fraction,
                             const double mz_bin_size);

      /// Statistics of how many neighbors were found per reference peptide
      struct NeighborStats
      {
        /** @name NeigborStats_members
         *  Mutually exclusive categories of how many neighbors were found per reference peptide
         */
        ///@{
        int unfindable_peptides = 0; ///< how many ref-peptides contain an 'X' (unknown amino acid) and thus cannot be searched for neighbors
        int findable_no_neighbors = 0; ///< how many peptides had no neighbors?
        int findable_one_neighbor = 0; ///< how many peptides had exactly one neighbor?
        int findable_multiple_neighbors = 0; ///< how many peptides had multiple neighbors?
        ///@} 
        
        /// Sum of all 4 categories
        int total() const
        {
          return unfindable_peptides + findable_no_neighbors + findable_one_neighbor + findable_multiple_neighbors;
        }
        /// Number of reference peptides that contain an 'X' (unknown amino acid), formatted as 'X (Y%)'
        String unfindable() const
        {
          return String(unfindable_peptides) + " (" + unfindable_peptides * 100 / total() + "%)";
        }

        /// Number of reference peptides that had no neighbors, formatted as 'X (Y%)'
        String noNB() const
        {
          return String(findable_no_neighbors) + " (" + findable_no_neighbors * 100 / total() + "%)";
        }
        /// Number of reference peptides that had exactly one neighbor, formatted as 'X (Y%)'
        String oneNB() const
        {
          return String(findable_one_neighbor) + " (" + findable_one_neighbor * 100 / total() + "%)";
        }
        /// Number of reference peptides that had multiple neighbors, formatted as 'X (Y%)'
        String multiNB() const
        {
          return String(findable_multiple_neighbors) + " (" + findable_multiple_neighbors * 100 / total() + "%)";
        }
      };

      /// after calling isNeighborPeptide() multiple times, this function returns the statistics of how many neighbors were found per reference peptide
      NeighborStats getNeighborStats() const;

    protected:
      /**
       * @brief Creates a map of masses to positions from the internal relevant peptides.
       * @return A map where the key is the mass and the value is a vector of positions.
       */
      std::map<double, std::vector<int>> createMassLookup_();
      
      /**
       * @brief Finds candidate positions based on a given mono-isotopic weight and mass tolerance.
       * @param mono_weight The mono-isotopic weight to find candidates for.
       * @param mass_tolerance The allowed tolerance for matching the mass.
       * @param mass_tolerance_pc_ppm Whether the mass tolerance is in ppm.
       * @return A pair of begin/end iterators into mass_position_map_ for the candidate positions
       */
      auto findCandidatePositions_(const double mono_weight, double mass_tolerance, const bool mass_tolerance_pc_ppm);


    private:
      const std::vector<AASequence>& digested_relevant_peptides_; ///< digested relevant peptides
      std::map<double, std::vector<int>> mass_position_map_; ///< map of masses to positions in digested_relevant_peptides_

      TheoreticalSpectrumGenerator spec_gen_; ///< for b/y ions with charge 1
      const Residue* x_residue_; ///< residue for unknown amino acid

      std::vector<int> neighbor_stats_; ///< how many neighbors per reference peptide searched using isNeighborPeptide()?

  }; // class NeighborSeq

} // namespace OpenMS
