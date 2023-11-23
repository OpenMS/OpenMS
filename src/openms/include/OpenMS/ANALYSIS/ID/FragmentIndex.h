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
      AASequence sequence; // containing the sequence including modification information
      UInt32 protein_idx;  // the Id from the Fasta file
      float precursor_mz;  // charge 1 precursor m/z
    };

    /// DefaultConstructor
    FragmentIndex();

    /// Default destructor
    ~FragmentIndex() override = default;

    // returns true if already built, false otherwise
    bool isBuild() const;

    enum class IonTypes
    {
      AX_IONS,
      BY_IONS,
      CZ_IONS
    };

    IonTypes getIonTypes() const;

    const std::vector<Peptide>& getPeptides() const;

#ifdef DEBUG_FRAGMENT_INDEX
    /**
     * @brief (only for debugging) Manually add a peptide to the DB. The peptide can be modified
     * IMPORTANT: add peptide after the other peptides are generated and before the FragmentIndex is build
     * @param peptide
     */
    void addSpecialPeptide(AASequence& peptide, Size source_idx);
#endif

    /** @brief Given a set of Fasta files, builds the Fragment Index datastructure (FID). First all fragments are sorted
     * by their own mass. Next they are placed in buckets. The min-fragment mass is stored for each bucket, whereupon
     * the fragments are sorted within the buckets by their originating precursor mass.
     *
     * @param fasta_entries
     */
    virtual void build(const std::vector<FASTAFile::FASTAEntry> & fasta_entries);

    /** @brief Delete fragment index. Sets is_build=false*/
    void clear();


    /** Return index range of all possible Peptides/Proteins, such that a vector can be created fitting that range (safe some memory)
     * @param precursor_mass The mono-charged precursor mass (M+H)
     * @param window Defines the lower and upper bound for the precusor mass. For closed search it only contains the tolerance. In case of open search
     *                  it contains both tolerance and open-search-window
     * @return a pair of indexes defining all possible peptides which the current peak could hit
     */
    std::pair<size_t, size_t> getPeptidesInPrecursorRange(float precursor_mass, std::pair<float, float> window);

    /**
     * A match between a single query peak and a database fragment
     */
    struct Hit
    {
      Hit(UInt32 peptide_idx, float fragment_mz) :
        peptide_idx(peptide_idx),
        fragment_mz(fragment_mz)
      {}
      UInt32 peptide_idx; // index in database
      float fragment_mz;
    };

    std::vector<Hit> query(Peak1D peak, std::pair<size_t, size_t> peptide_idx_range, uint16_t peak_charge);

protected:
    bool is_build_{false};              ///< true, if the database has been populated with fragments

    void updateMembers_() override;

     /**@brief Generates all peptides from given fasta entries. If Bottom-up is set to false
     * skips digestion. If set to true the Digestion enzyme can be set in the parameters.
     * Additionally introduces fixed and variable modifications for restrictive PSM search.
     *
     * @param fasta_entries
     */
    void generate_peptides(const std::vector<FASTAFile::FASTAEntry>& fasta_entries);

    std::vector<Peptide> fi_peptides_;   ///< vector of all (digested) peptides

    float fragment_min_mz_;  ///< smallest fragment mz
    float fragment_max_mz_;  ///< largest fragment mz    
    size_t bucketsize_;       ///< number of fragments per outer node
    std::vector<float> bucket_min_mz_;  ///< vector of the smalles fragment mz of each bucket
    float precursor_mz_tolerance_;
    bool precursor_mz_tolerance_unit_ppm_{true};
    float fragment_mz_tolerance_;
    bool fragment_mz_tolerance_unit_ppm_{true};    
private:
  
    /**@brief One entry in the fragment index
     */
    struct Fragment
    {
      UInt32 peptide_idx; // TODO: check size in sage implementation (32bit vs 64bit)
      float fragment_mz;
    };

    //since we work with TheoreticalSpectrumGenerator, we must transfer some of those member variables
    bool add_b_ions_;
    bool add_y_ions_;
    bool add_a_ions_;
    bool add_c_ions_;
    bool add_x_ions_;
    bool add_z_ions_;

    // SpectrumGenerator independend member variables
    std::string digestion_enzyme_;
    IonTypes ion_types_{FragmentIndex::IonTypes::BY_IONS};

    size_t missed_cleavages_; ///< number of missed cleavages
    float peptide_min_mass_;
    float peptide_max_mass_;
    size_t peptide_min_length_;
    size_t peptide_max_length_;
  
    StringList modifications_fixed_;    ///< Modification that are one all peptides
    StringList modifications_variable_; ///< Variable Modification -> all possible comibnations are created
    size_t max_variable_mods_per_peptide_;

    std::vector<Fragment> fi_fragments_; ///< vector of all theoretical fragments (b- and y- ions)
  };

}
