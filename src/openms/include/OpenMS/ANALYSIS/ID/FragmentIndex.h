// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/DATASTRUCTURES/MultiFragment.h>

#include <vector>
#include <functional>


namespace OpenMS
{
  /** @brief Generates from a set of Fasta files a 2D-datastructure which stores all theoretical masses of all
   * b and y ions from all peptides generated from the Fasta file. The datastructure is build such that on one axis
   * the fragments are sorted by their own mass and the axis by the mass of their precursor/protein.
   * The FI has two options: Bottom-up and Top Down. In later digestion is skiped and the fragments have a direct
   * reference to the mass of the proteins instead of digested peptides.
   */
  class OPENMS_DLLAPI FragmentIndex : public DefaultParamHandler
  {


  public:

    /** @brief Peptide with all important infos needed for the FI-structure
     */
    struct Peptide
    {
      bool expType;        //0 for top-down or 1 for bottom-up
      AASequence sequence; // containing the sequence including modification information
      Size protein_idx;    // the Id from the Fasta file
      double mass;         // precursor-mass

    };


    /// DefaultConstructor
    FragmentIndex();

    /// Default destructor
    ~FragmentIndex() override = default;

    /// assignment operator
    virtual FragmentIndex& operator=(const FragmentIndex& fi);

    //getter
    std::vector<AASequence> getFiPeptidesSequences() const;
    bool isBuild() const;
    bool containsFragmentation(const std::string& fragmentation) const;

    const std::vector<Peptide>& getFiPeptides() const;

    /**
     * @brief (only for debugging) Manually add a peptide to the DB. The peptide can be modified
     * IMPORTANT: add peptide after the other peptides are generated and before the FragmentIndex is build
     * @param peptide
     */
    void addSpecialPeptide(AASequence& peptide, Size source_idx);

    void clear();


    /**@brief Generates all peptides from given fasta entries. If Bottom-up is set to false
     * skips digestion. If set to true the Digestion enzyme can be set in the parameters.
     * Additionally introduces fixed and variable modifications for restrictive PSM search.
     *
     * @param fasta_entries
     */
    void generate_peptides(const std::vector<FASTAFile::FASTAEntry> & fasta_entries);

    /** @brief Given a set of Fasta files, builds the Fragment Index datastructure (FID). First all fragments are sorted
     * by their own mass. Next they are placed in buckets. The min-fragment mass is stored for each bucket, whereupon
     * the fragments are sorted within the buckets by their originating precursor mass.
     *
     * @param fasta_entries
     */
    virtual void build(const std::vector<FASTAFile::FASTAEntry> & fasta_entries);


    /** General function for binary-search given a upper and lower bound. Is needed to query for the
     * Precursor mass, as well for the fragment mass
     *
     * @tparam S The class that is stored inside the vector, which we want to search
     * @tparam B The class that we want to extract from each class S, to make a comparison
     * @param slice The input vector
     * @param low Lower bound value
     * @param high Higher bound value
     * @param access A lambda function that defines how B can be accessed from S
     * @param include Whether the found lower bound element should be included or not (important when we search the fragment masses
     *                  where not the actual mass is stored but the min-mass
     * @return a pair of indexes containing the lower bound and upper bound
     */
    template<class S, class B> static std::pair<size_t, size_t> binary_search_slice(const std::vector<S>& slice, B low, B high, B (*access)(S), bool include);

    static std::pair<size_t, size_t > binary_search_slice_double(const std::vector<double>& slice, double low, double high, bool include);
    static std::pair<size_t, size_t > binary_search_slice_mf(const std::vector<OpenMS::MultiFragment>& slice,
                                                         size_t low, size_t high,
                                                         size_t (*access) (OpenMS::MultiFragment),
                                                         bool include);

    /**
     * A match between query peak and database peak
     */
    struct Hit
    {
      size_t peptide_idx;
      double fragment_mz;
    };

    /** Return index range of all possible Peptides/Proteins, such that a vector can be created fitting that range (safe some memory)
     *
     * @param precursor_mass The precursor mass we search for
     * @param window Defines the lower and upper bound for the precusor mass. For closed search it only contains the tolerance. In case of open search
     *                  it contains both tolerance and open-search-window
     * @return a pair of indexes defining all possible peptides which the current peak could hit
     */
    std::pair<size_t, size_t > getPeptideRange(double precursor_mass, std::pair<double, double> window);

    std::vector<Hit> query(Peak1D peak, std::pair<size_t, size_t> peptide_idx_range, uint16_t peak_charge);



  protected:

    void updateMembers_() override;


    /**@brief The struct the Fragments of the FID are safed in .
     *
     */
    struct Fragment
    {
      Size peptide_idx;
      double fragment_mz;
      //double precursor_mz;  // this one is not needed, bc the peptide vector is sorted, so we simply look for the peptide_idx

    };

  protected:

    //since we work with TheoreticalSpectrumGenerator, we must transfer some of those member variables
    bool add_b_ions_;
    bool add_y_ions_;
    bool add_a_ions_;
    bool add_c_ions_;
    bool add_x_ions_;
    bool add_z_ions_;

    // SpectrumGenerator independend member variables
    bool exp_type_;   ///< Bottom-up or Top-Down
    std::string digestion_enzyme_;
    size_t missed_cleavages_; ///< number of missed cleavages
    double peptide_min_mass_;
    double peptide_max_mass_;
    size_t peptide_min_length_;
    size_t peptide_max_length_;
    double fragment_min_mz_;
    double fragment_max_mz_;
    size_t bucketsize_; ///< number of fragments per outer node
    double precursor_mz_tolerance_;
    std::string precursor_mz_tolerance_unit_;
    double fragment_mz_tolerance_;
    std::string fragment_mz_tolerance_unit_;
    StringList modifications_fixed_; ///< Modification that are one all peptides
    StringList modifications_variable_; ///< Variable Modification -> all possible comibnations are created
    size_t max_variable_mods_per_peptide_;
    bool is_build_;               ///< true, if the database is built
    std::vector<Peptide> fi_peptides_;   ///< vector of all (digested) peptides
    std::vector<Fragment> fi_fragments_; ///< vector of all theoretical fragments (b- and y- ions)
    std::vector<double> bucket_min_mz_;  ///< vector of the smalles fragment mz of each bucket

  };

}