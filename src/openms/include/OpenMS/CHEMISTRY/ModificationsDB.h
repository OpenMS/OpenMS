// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_MODIFICATIONSDB_H
#define OPENMS_CHEMISTRY_MODIFICATIONSDB_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <set>

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
      databases. This can be done by providing a path when initializing
      ModificationsDB.
  */
  class OPENMS_DLLAPI ModificationsDB
  {
public:

    /// Returns a pointer to the modifications DB (singleton)
    inline static ModificationsDB* getInstance(OpenMS::String unimod_file = "CHEMISTRY/unimod.xml", OpenMS::String psimod_file = "CHEMISTRY/PSI-MOD.obo", OpenMS::String xlmod_file = "CHEMISTRY/XLMOD.obo")
    {
      static ModificationsDB* db_ = nullptr;
      if (db_ == nullptr)
      {
        db_ = new ModificationsDB(unimod_file, psimod_file, xlmod_file);
      }
      return db_;
    }

    /// Check whether ModificationsDB was instantiated before
    static bool isInstantiated();

    friend class CrossLinksDB;

    /// Returns the number of modifications read from the unimod.xml file
    Size getNumberOfModifications() const;

    /**
       @brief Returns the modification with the given index

       @throw Exception::IndexOverflow if the index is too large
    */
    const ResidueModification& getModification(Size index) const;

    /**
       @brief Collects all modifications which have the given name as synonym

       If @p residue is set, only modifications with matching residue of origin are considered.
       If @p term_spec is set, only modifications with matching term specificity are considered.
       The resulting set of modifications may be empty if no modification exists that fulfills the criteria.

       @throw Exception::ElementNotFound if no modification named @p mod_name exists
       @throw Exception::InvalidValue if no residue named @p residue exists
    */
    void searchModifications(std::set<const ResidueModification*>& mods, const String& mod_name, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY) const;

    /**
       @brief Returns the modification with the given name

       If @p residue is set, only modifications with matching residue of origin are considered.
       If @p term_spec is set, only modifications with matching term specificity are considered.

       If more than one matching modification is found, the first one is returned with a warning.

       @throw Exception::ElementNotFound if no modification named @p mod_name exists (via searchModifications())
       @throw Exception::InvalidValue if no matching modification exists
    */
    const ResidueModification& getModification(const String& mod_name, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY) const;

    /// Returns true if the modification exists
    bool has(String modification) const;

    /**
       @brief Add a new modification to ModificationsDB.
       @throw Exception::InvalidValue if modification already exists (based on its fullID)
    */
    void addModification(ResidueModification * new_mod);

    /**
       @brief Returns the index of the modification in the mods_ vector; a unique name must be given

       @throw Exception::ElementNotFound if not exactly one matching modification was found
       @throw Exception::InvalidValue if no matching residue or modification were found
    */
    Size findModificationIndex(const String& mod_name) const;

    /**
       @brief Collects all modifications with matching delta mass

       If @p residue is set, only modifications with matching residue of origin are considered.
       If @p term_spec is set, only modifications with matching term specificity are considered.
    */
    void searchModificationsByDiffMonoMass(std::vector<String>& mods, double mass, double max_error, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY);

    /** @brief Returns the best matching modification for the given mass and residue

        Query the modifications DB to get the best matching modification with
        the given mass at the given residue (NULL pointer means no result,
        maybe the maximal error tolerance needs to be increased). Possible
        input for CAM modification would be a mass of 160 and a residue of
        "C".

        @note If there are multiple possible matches with equal masses, it
        will choose the _first_ match which defaults to the first matching
        UniMod entry.

        @param residue The residue at which the modifications occurs
        @param mass The monoisotopic mass of the residue including the mass of the modification
        @param max_error The maximal mass error in the modification search

        @return A pointer to the best matching modification (or NULL if none was found)
    */
    const ResidueModification* getBestModificationByMonoMass(double mass, double max_error, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY);

    /** @brief Returns the best matching modification for the given delta mass and residue

        Query the modifications DB to get the best matching modification with
        the given delta mass at the given residue (NULL pointer means no result,
        maybe the maximal error tolerance needs to be increased). Possible
        input for CAM modification would be a delta mass of 57 and a residue
        of "C".

        @note If there are multiple possible matches with equal masses, it
        will choose the _first_ match which defaults to the first matching
        UniMod entry.

        @param residue The residue at which the modifications occurs
        @param mass The monoisotopic mass of the residue including the mass of the modification
        @param max_error The maximal mass error in the modification search

        @return A pointer to the best matching modification (or NULL if none was found)

    */
    const ResidueModification* getBestModificationByDiffMonoMass(double mass, double max_error, const String& residue = "", ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY);

    /// Collects all modifications that can be used for identification searches
    void getAllSearchModifications(std::vector<String>& modifications) const;

protected:

    /// Stores whether ModificationsDB was instantiated before
    static bool is_instantiated_;

    /// Stores the modifications
    std::vector<ResidueModification*> mods_;

    /// Stores the mappings of (unique) names to the modifications
    Map<String, std::set<const ResidueModification*> > modification_names_;

    /// Helper function to check if a residue matches the origin for a modification
    bool residuesMatch_(const String& residue, char origin) const;

private:

    /** @name Constructors and Destructors

        @param unimod_file Path to the Unimod XML file
        @param psimod_file Path to the PSI-MOD OBO file
        @param xlmod_file Path to the XLMOD OBO file

     */
    //@{
    /// Default constructor
    ModificationsDB(OpenMS::String unimod_file = "CHEMISTRY/unimod.xml", OpenMS::String psimod_file = "CHEMISTRY/PSI-MOD.obo", OpenMS::String xlmod_file = "CHEMISTRY/XLMOD.obo");

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
       @brief Adds modifications from a given file in OBO format

       @throw Exception::ParseError if the file cannot be parsed correctly
    */
    void readFromOBOFile(const String& filename);

    /// Adds modifications from a given file in Unimod XML format
    void readFromUnimodXMLFile(const String& filename);
    
  };
}
#endif
