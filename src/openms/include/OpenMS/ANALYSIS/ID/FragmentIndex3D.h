// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGeneratorNode.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/MultiFragment.h>
#include <OpenMS/DATASTRUCTURES/MultiPeak.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <functional>
#include <utility>
#include <vector>


namespace OpenMS
{
  /** @brief Generates from a set of Fasta files a 3D-datastructure which stores all theoretical masses of all
   * b and y ions from all peptides generated from the Fasta file. The datastructure is build such that on one axis
   * the fragments are sorted by their own mass and the axis by the mass of their precursor/protein.
   * In addition a third axis is introduced where the Fragments are sorted based on neighbor infomation.
   * The FI has two options: Bottom-up and Top Down. In later digestion is skiped and the fragments have a direct
   * reference to the mass of the proteins instead of digested peptides.
   */
  class OPENMS_DLLAPI FragmentIndex3D : public FragmentIndex
  {


  public:


    /**
     * Builds the Index Database in a multi-dimensional/multi-leveled tree structure
     * @param fasta_entries
     */
    void build(const std::vector<FASTAFile::FASTAEntry> &fasta_entries) override;

    /**
     * Takes a MultiPeak as input and recursively searches for a hit in the DB
     * @param hits output
     * @param peak the query Multipeak
     * @param peptide_idx_range The range of all possible Peptides
     * @param window The window in which we want to search, enabeling finding of modified peptides
     */
    void query(std::vector<Hit>& hits, const MultiPeak& peak, std::pair<size_t, size_t> peptide_idx_range, std::pair<double, double> window);


    /// DefaultConstructor
    FragmentIndex3D() ;

    /// Default destructor
    //FragmentIndex;

  protected:

    void updateMembers_() override;

  private:

    /**
     * @brief Check if potential hit is in the range of [query-tolerance + window, query + tolerance + window]
     * @param hit entry in the DB which could be a potential hit
     * @param query the current query
     * @param tolerance
     * @param window
     * @return
     */
    static bool inRange(double hit, double query, double tolerance, std::pair<double, double> window);
    /**
     * @brief compares two vectors and checks if each single element of the hit vector is in range of its corresponding query counterpart
     * @param hit
     * @param query
     * @param tolerance
     * @return
     */
    static bool inRangeFollowUpPeaks(std::vector<double> hit, std::vector<double> query, double tolerance);

    uint16_t depth_; // The depth of the database (e.q. Depth 3. We include the next three peaks on the right. The database is then (3+2) Dimensional)
    std::vector<MultiFragment> fi_fragments_;
    std::vector<std::vector<double>> follow_up_peaks_buckets_min_mz;

    void recursiveQuery(std::vector<Hit>& hits,
                        const MultiPeak& peak,
                        std::pair<size_t, size_t> peptide_idx_range,
                        std::pair<double, double> window,  // The window for the Fragment mass ONLY!!!
                        size_t recursion_step,
                        size_t current_slice,
                        double fragment_tolerance);



  };

}