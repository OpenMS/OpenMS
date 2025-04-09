// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <set>
#include <memory>  // unique_ptr
#include <unordered_map>

namespace OpenMS
{
  // forward declarations
  class ResidueModification;
  class Residue;

  /** @ingroup Chemistry

      @brief database which holds all residue modifications from UniMod

      This singleton class serves as a storage of the available modifications
      represented by UniMod (www.unimod.org). The modifications are identified
      by their name and possibly other IDs from UniMod or the PSI-MOD ontology.
      Modifications can have different specificities, e.g. they can occur only
      at the termini, anywhere or only at specific amino acids.

      The modifications are defined in share/OpenMS/CHEMISTRY/unimod.xml and
      in share/OpenMS/CHEMISTRY/PSI-MOD.obo. The unimod file can be directly
      downloaded from unimod.org and replaced if the modifications change.

      To add a new modification, not contained in UniMod, one should follow
      the way described at the unimod.org website and download the file then
      from unimod.org. The same can be done to add support for the modifications
      to search engines, e.g. Mascot.

      In some scenarios, it might be useful to define different modification
      databases. This can be done by providing a path through
      initializeModificationsDB(), however it is important that this is done
      *before* the first call to getInstance().
  */
  class OPENMS_DLLAPI ModificationsDB
  {
public:

    /// Returns a pointer to the modifications DB (singleton)
    static ModificationsDB* getInstance();

    /// Initializes the modification DB with non-default modification files (can only be done once)
    static ModificationsDB* initializeModificationsDB(OpenMS::String unimod_file = "CHEMISTRY/unimod.xml", OpenMS::String custommod_file = "CHEMISTRY/custom_mods.xml", OpenMS::String psimod_file = "CHEMISTRY/PSI-MOD.obo", OpenMS::String xlmod_file = "CHEMISTRY/XLMOD.obo");

    /// Check whether ModificationsDB was instantiated before
    static bool isInstantiated();

    friend class CrossLinksDB;
    // for access to addNewModification_ (without checking presence)
    friend class Residue;
    friend class AASequence;

    /// Returns the number of modifications read from the unimod.xml file
    Size getNumberOfModifications() const;

    /**
       @brief Returns the modification with the given index.
       note: out-of-bounds check is only performed in debug mode.
    */
    const ResidueModification* getModification(Size index) const;

    /**
       @brief Collects all modifications which have the given name as synonym
       @todo use set as return value. Would be more efficient in pyopenms

       If @p residue is set, only modifications with matching residue of origin are considered.
       If @p term_spec is set, only modifications with matching term specificity are considered.
        The resulting set of modifications will be empty if no modification exists that fulfills the criteria.
    */
    void searchModifications(std::set<const ResidueModification*>& mods,
                             const String& mod_name,
                             const String& residue = "",
                             ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY) const;

   /**
      @brief Returns a pointer to an exact match of the given modification if present in the DB.

      This should be used if e.g. only a stack copy of the modification is available but you need
      a pointer to the modification in the database.

      @return The matching modification given the constraints. Returns nullptr
      if no modification exists that is an exact match accoring to the equals operator.
   */
   const ResidueModification* searchModification(const ResidueModification& mod_in) const;

    /**
       @brief Returns the modification which has the given name as synonym (fast version)

       Unlike searchModifications(), only returns the one occurrence of the
       modification (the last occurrence). It is therefore required to check @p
       multiple_matches to ensure that only a single modification was found.

       If @p residue is set, only modifications with matching residue of origin are considered.
       If @p term_spec is set, only modifications with matching term specificity are considered.

       @return The matching modification given the constraints. Returns nullptr
       if no modification exists that fulfills the criteria. If multiple
       modifications are found, the @p multiple_matches flag will be set.
    */
    const ResidueModification* searchModificationsFast(const String& mod_name,
                                                       bool& multiple_matches,
                                                       const String& residue = "",
                                                       ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY) const;

    /**
       @brief Returns the modification with the given name

       If @p residue is set, only modifications with matching residue of origin are considered.
       If @p term_spec is set, only modifications with matching term specificity are considered.

       If more than one matching modification is found, the first one is returned with a warning.

        @note Will never return a null pointer, instead will throw an exceptions.
       
        @throw Exception::ElementNotFound if no modification named @p mod_name exists (via searchModifications())
       @throw Exception::InvalidValue if no matching modification exists
    */
    const ResidueModification* getModification(const String& mod_name, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY) const;

    /// Returns true if the modification exists
    bool has(const String& modification) const;

    /**
       @brief Add a new modification to ModificationsDB.
       If the modification already exists (based on its fullID) it is not added.
       @return a pointer to the modification in the ModificationDB (which can differ from input if mod was already present).

       @param new_mod Owning pointer, which transfers ownership to ModificationsDB (mod might get deleted if already present!)
    */
    const ResidueModification* addModification(std::unique_ptr<ResidueModification> new_mod);

    /**
       @brief Add a new modification to ModificationsDB.
       If the modification already exists (based on its fullID) it is not added. A copy will be made on the heap and added to the ModificationsDB otherwise.
       @return a pointer to the modification in the ModificationDB (which can differ from input if mod was already present).

       @param new_mod The new modification object. A copy will be made on the heap and added to the ModificationsDB if not already present.
    */
    const ResidueModification* addModification(const ResidueModification& new_mod);

    /**
       @brief Returns the index of the modification in the mods_ vector; a unique name must be given

       return numeric_limits<Size>::max() if not exactly one matching modification was found
       or no matching residue or modification were found

       @throw Exception::ElementNotFound if not exactly one matching modification was found. 
    */
    Size findModificationIndex(const String& mod_name) const;

    /**
       @brief Collects all modifications with delta mass inside a tolerance window
       @warning This function adds the results in the order of appearance in the DB, not considering proximity in mass. Use searchModificationsByDiffMonoMassSorted for this.

       If @p residue is set, only modifications with matching residue of origin are considered.
       If @p term_spec is set, only modifications with matching term specificity are considered.
    */
    void searchModificationsByDiffMonoMass(std::vector<String>& mods, double mass, double max_error, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
    void searchModificationsByDiffMonoMass(std::vector<const ResidueModification*>& mods, double mass, double max_error, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY);

    /**
     @brief Collects all modifications with delta mass inside a tolerance window and adds them sorted
     by mass difference

     If @p residue is set, only modifications with matching residue of origin are considered.
     If @p term_spec is set, only modifications with matching term specificity are considered.
    */
    void searchModificationsByDiffMonoMassSorted(std::vector<String>& mods, double mass, double max_error, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
    void searchModificationsByDiffMonoMassSorted(std::vector<const ResidueModification*>& mods, double mass, double max_error, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY);


    /** @brief Returns the best matching modification for the given delta mass and residue

        Query the modifications DB to get the best matching modification with
        the given delta mass at the given residue (NULL pointer means no result,
        maybe the maximal error tolerance needs to be increased). Possible
        input for CAM modification would be a delta mass of 57 and a residue
        of "C".

        @note If there are multiple possible matches with equal masses, it
        will choose the _first_ match which defaults to the first matching
        UniMod entry.

        @param mass The monoisotopic mass of the residue including the mass of the modification
        @param max_error The maximal mass error in the modification search
        @param residue The residue at which the modifications occurs
        @param term_spec Only modifications with matching term specificity are considered.

        @return A pointer to the best matching modification (or NULL if none was found)

    */
    const ResidueModification* getBestModificationByDiffMonoMass(double mass, double max_error, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY);

    /// Collects all modifications that can be used for identification searches
    void getAllSearchModifications(std::vector<String>& modifications) const;

    /// Writes tab separated entries: FullId,FullName,Origin,AA,TerminusSpecificity,DiffMonoMass (including header) to TSV file
    void writeTSV(const String& filename);

  protected:

    /// Stores whether ModificationsDB was instantiated before
    static bool is_instantiated_;

    /// Stores the modifications
    std::vector<ResidueModification*> mods_;

    /// Stores the mappings of (unique) names to the modifications
    std::unordered_map<String, std::set<const ResidueModification*> > modification_names_;

    /** @brief Helper function to check if a residue matches the origin for a modification
     *
     * Special cases are handled as follows:
     *   * if the origin of the modification is not "X" (everything), then the
     *     residue either needs to match the origin exactly or it must be one of "X", ".", or "?"
     *   * if the origin of the modification is "X" (can match any amino acid),
     *     then any residue should match -- except if the modification is
     *     user-defined and maps to an unknown amino acid (indicated by "X")
     *
     * Underlying logic to determine whether a given residue matches the
     * modification: if the modification does not have origin of "X"
     * (everything) then it is sufficient to check that the residue matches the origin
     *
    */
    bool residuesMatch_(const char residue, const ResidueModification* curr_mod) const;

private:

    /** @name Constructors and Destructors

        @param unimod_file Path to the Unimod XML file
        @param psimod_file Path to the PSI-MOD OBO file
        @param xlmod_file Path to the XLMOD OBO file

     */
    //@{
    explicit ModificationsDB(const OpenMS::String& unimod_file = "CHEMISTRY/unimod.xml", const OpenMS::String& custommod_file = "CHEMISTRY/custom_mods.xml", const OpenMS::String& psimod_file = "CHEMISTRY/PSI-MOD.obo", const OpenMS::String& xlmod_file = "CHEMISTRY/XLMOD.obo");

    /// Copy constructor
    ModificationsDB(const ModificationsDB& residue_db);

    /// Destructor
    virtual ~ModificationsDB();
    //@}

    /** @name Assignment
     */
    //@{
    /// Assignment operator
    ModificationsDB & operator=(const ModificationsDB& aa);
    //@}

    /**
       @brief Add a new modification to ModificationsDB without checking if it was inside already.

       @param new_mod A copy will be made on the heap and added to the modification if not already present.
    */
    const ResidueModification* addNewModification_(const ResidueModification& new_mod);

    /**
       @brief Adds modifications from a given file in OBO format

       @throw Exception::ParseError if the file cannot be parsed correctly
    */
    void readFromOBOFile(const String& filename);

    /// Adds modifications from a given file in Unimod XML format
    void readFromUnimodXMLFile(const String& filename);
  };
}
