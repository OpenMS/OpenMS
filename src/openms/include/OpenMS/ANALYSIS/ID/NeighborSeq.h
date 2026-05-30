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
      @brief Subset-neighbor peptide search: find peptides from a wider
             pool (typically a FASTA digest) that are spectral neighbors
             of a smaller "relevant" peptide set, useful when only part
             of a complex sample is of interest.

      Two peptides are considered neighbors when their precursor masses
      are within tolerance and their theoretical b/y fragment spectra
      share enough peaks. The class is constructed once with the
      relevant peptides, then queried with @ref isNeighborPeptide for
      each candidate. After the queries are done, @ref getNeighborStats
      summarises how many relevant peptides had zero, one, or multiple
      neighbors.

      Background: Cormen et al., @e J. Proteome Research 2021,
      @c 10.1021/acs.jproteome.1c00483
      (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8489664/).

      @ingroup Analysis_ID
    */
    class OPENMS_DLLAPI NeighborSeq
    {

    public:
      /**
        @brief Construct from a vector of "relevant" digested peptides.

        Builds an internal mass index (relevant peptides containing the
        unknown-amino-acid residue @c 'X' are skipped and a count is
        logged via OPENMS_LOG_WARN) and configures the internal
        theoretical-spectrum generator for charge-1 b/y ions including
        the @c b1 prefix ion.

        @note The class stores a @c const-reference to the moved-in
              vector via its internal member (not a copy). The vector
              passed to the constructor must therefore outlive every
              call on this instance.

        @param[in] digested_relevant_peptides Digested peptides to use
                                              as the "relevant" reference
                                              set.
      */
      NeighborSeq(std::vector<AASequence>&& digested_relevant_peptides);

      /**
       * @brief Generates a theoretical spectrum for a given peptide sequence with b/y ions at charge 1.
       * 
       * Includes all b and y ions with charge 1 (even the prefix ions, e.g. b1), but no internal ions.
       * 
       * @param[in] peptide_sequence The peptide sequence for which to generate the spectrum.
       * @return The generated theoretical spectrum.
       */
      MSSpectrum generateSpectrum(const AASequence& peptide_sequence);

      /**
        @brief Whether two spectra share enough peaks (in @p mz_bin_size m/z bins) to be considered neighbors.

        Computes the Dice-style fraction
        @c 2 @c * @c shared_peaks @c / @c (spec1.size() @c + @c spec2.size())
        and returns @c true when it is @b strictly greater than
        @p min_shared_ion_fraction. All peaks of both spectra are
        considered.

        @param[in] spec1                   First theoretical spectrum.
        @param[in] spec2                   Second theoretical spectrum.
        @param[in] min_shared_ion_fraction Minimum required Dice-style
                                           shared-peak fraction (in
                                           @c [0, @c 1]).
        @param[in] mz_bin_size             Bin size for m/z comparison
                                           (typical values: @c 0.05 Th
                                           for high-resolution data,
                                           @c 1.0005079 Th for
                                           low-resolution data).
        @return @c true when the shared fraction is strictly greater
                than @p min_shared_ion_fraction.
      */
      static bool isNeighborSpectrum(const MSSpectrum& spec1, const MSSpectrum& spec2, const double min_shared_ion_fraction, const double mz_bin_size);
      /**
       * @brief Compute the number of shared ions between two spectra 
       *
       * All peaks are considered. Use generateSpectrum() to generate theoretical spectra with b/y ions.
       * 
       * @param[in] spec1 The first theoretical spectrum.
       * @param[in] spec2 The second theoretical spectrum.
       * @param[in] mz_bin_size Bin size for the m/z values, which determines if two peaks are considered to be the same.
       * @return The number of shared ions
       */
      static int computeSharedIonCount(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& mz_bin_size);

      /**
        @brief Whether @p neighbor_candidate is a spectral neighbor of any of the relevant peptides.

        Looks up the relevant peptides whose precursor mass is within
        @p mass_tolerance_pc of @p neighbor_candidate's mono-isotopic
        mass, generates b/y spectra for each candidate match plus
        @p neighbor_candidate, and compares them with
        @ref isNeighborSpectrum. Returns @c true as soon as @b any
        relevant peptide qualifies, but continues iterating so that the
        internal per-relevant-peptide neighbor counters are updated for
        every match. Call @ref getNeighborStats once all candidates
        have been queried.

        @param[in] neighbor_candidate     Candidate peptide (typically
                                          from a digested FASTA).
        @param[in] mass_tolerance_pc      Maximum precursor mass
                                          difference between
                                          @p neighbor_candidate and a
                                          relevant peptide, expressed
                                          in Da or ppm per
                                          @p mass_tolerance_pc_ppm.
        @param[in] mass_tolerance_pc_ppm  @c true to interpret
                                          @p mass_tolerance_pc as ppm
                                          (converted internally via
                                          @c Math::ppmToMass);
                                          @c false to interpret it as
                                          Da.
        @param[in] min_shared_ion_fraction Threshold passed straight to
                                           @ref isNeighborSpectrum.
        @param[in] mz_bin_size            m/z bin size passed straight
                                          to @ref isNeighborSpectrum
                                          (typical values: @c 0.05 Th
                                          for high-res, @c 1.0005079 Th
                                          for low-res).
        @return @c true when at least one relevant peptide is a neighbor.
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
        
        /// Sum of all four categories (i.e. the number of relevant peptides registered at construction time).
        int total() const
        {
          return unfindable_peptides + findable_no_neighbors + findable_one_neighbor + findable_multiple_neighbors;
        }

        /**
          @brief @ref unfindable_peptides formatted as @c "X (Y%)".

          @warning Triggers integer division by zero when @ref total is @c 0
                   (the four counters and the formatter share an integer denominator).
        */
        String unfindable() const
        {
          return String(unfindable_peptides) + " (" + unfindable_peptides * 100 / total() + "%)";
        }

        /// @ref findable_no_neighbors formatted as @c "X (Y%)"; see @ref unfindable for the divide-by-zero caveat.
        String noNB() const
        {
          return String(findable_no_neighbors) + " (" + findable_no_neighbors * 100 / total() + "%)";
        }

        /// @ref findable_one_neighbor formatted as @c "X (Y%)"; see @ref unfindable for the divide-by-zero caveat.
        String oneNB() const
        {
          return String(findable_one_neighbor) + " (" + findable_one_neighbor * 100 / total() + "%)";
        }

        /// @ref findable_multiple_neighbors formatted as @c "X (Y%)"; see @ref unfindable for the divide-by-zero caveat.
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
       * @param[in] mono_weight The mono-isotopic weight to find candidates for.
       * @param[in] mass_tolerance The allowed tolerance for matching the mass.
       * @param[in] mass_tolerance_pc_ppm Whether the mass tolerance is in ppm.
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
