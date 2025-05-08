// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Jang Jang Jin$
// --------------------------------------------------------------------------

#pragma once

#include <unordered_map>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Macros.h> // for OPENMS_PRECONDITION

#include <map>
#include <set>
#include <array>

namespace OpenMS
{
  // forward declarations
  class ResidueModification;
  class Residue;

  /** @ingroup Chemistry

      @brief OpenMS stores a central database of all residues in the ResidueDB.
      All (unmodified) residues are added to the database on construction.
      Modified residues get created and added if getModifiedResidue is called.
  */
  class OPENMS_DLLAPI ResidueDB
  {
public:
    /// singleton
    static ResidueDB* getInstance();

    /** @name Constructors and Destructors
    */
    //@{
    /// destructor
    virtual ~ResidueDB();
    //@}

    /** @name Accessors
    */
    //@{
    /// returns the number of residues stored
    Size getNumberOfResidues() const;

    /// returns the number of modified residues stored
    Size getNumberOfModifiedResidues() const;

    /// returns a pointer to the residue with name, 3 letter code or 1 letter code name
    const Residue* getResidue(const String& name) const;

    /// returns a pointer to the residue with 1 letter code name
    const Residue* getResidue(const unsigned char& one_letter_code) const;

    /**
       @brief Returns a pointer to a modified residue given a modification name

       The "base" residue is looked up in ModificationsDB using the modification name.
       The modified residue is added to the database if it doesn't exist yet.
    */
    const Residue* getModifiedResidue(const String& name);

    /**
       @brief Returns a pointer to a modified residue given a residue and a modification name

       The modified residue is added to the database if it doesn't exist yet.

       @throw Exception::IllegalArgument if the residue was not found
       @throw Exception::InvalidValue if no matching modification was found (via ModificationsDB::getModification)
    */
    const Residue* getModifiedResidue(const Residue* residue, const String& name);

    /**
       @brief Returns a pointer to a modified residue given a residue and a pointer to a modification from the ModificationsDB

       The modified residue is added to the database if it doesn't exist yet. The origin of the modification has to match the residue and
       the term has to be ResidueModification::Anywhere.

       @throw Exception::IllegalArgument if the residue was not found
    */
    const Residue* getModifiedResidue(const Residue* residue, const ResidueModification* mod);

    /**
       @brief returns a set of all residues stored in this residue db

       Following sets are available:
       All - all residues
       Natural20 - default 20 naturally occurring residues
       Natural19WithoutI - default natural amino acids, excluding isoleucine (isobaric to leucine)
       Natural19WithoutL - default natural amino acids, excluding leucine (isobaric to isoleucine)
       Natural19J - default natural amino acids,  (isobaric leucine/isoleucine are marked by 'J')
       AmbiguousWithoutX - all amino acids, including ambiguous ones: B (asparagine or aspartate), Z (glutamine or glutamate), J (isoleucine or leucine)
       Ambiguous - all amino acids including all ambiguous ones (X can be every other amino acid)
       AllNatural - naturally occurring residues, including selenocysteine (U)

       returns an empty set if the specified residue set is not defined
    */
    const std::set<const Residue*> getResidues(const String& residue_set = "All") const;

    /// returns all residue sets that are registered which this instance
    const std::set<String> getResidueSets() const;

    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if the db contains a residue with the given name
    bool hasResidue(const String& name) const;

    /// returns true if the db contains the residue of the given pointer
    bool hasResidue(const Residue* residue) const;
    //@}

protected:
    /// initializes all residues by building
    void initResidues_();

    /** @name Private Constructors
    */
    //@{
    /// default constructor
    ResidueDB();

    ///copy constructor
    ResidueDB(const ResidueDB& residue_db);
    //@}

    /** @name Assignment
    */
    //@{
    /// assignment operator
    ResidueDB& operator=(const ResidueDB& aa) = delete;
    //@}

   // construct all residues 
    void buildResidues_();
    
    /// creates and adds residues to a lookup table including the residue set
    void insertResidueAndAssociateWithResidueSet_(Residue* residue, const std::vector<String>& residue_sets);

    /// add residue and add names to lookup
    void addResidue_(Residue* residue);

    /// adds names of single residue to the index
    void addResidueNames_(const Residue*);

    /// adds names of single modified residue to the index
    void addModifiedResidueNames_(const Residue*);
    
    std::map<String, std::map<String, const Residue*> > residue_mod_names_;

    /// all (unmodified) residues
    std::set<const Residue*> const_residues_;

    /// all modified residues
    std::set<const Residue*> const_modified_residues_;

    std::set<String> residue_sets_;

    /// lookup from name to residue
    std::unordered_map<String, const Residue*> residue_names_;

    /// fast lookup table for residues  
    std::array<const Residue*, 256> residue_by_one_letter_code_ = {{nullptr}};

    std::map<String, std::set<const Residue*> > residues_by_set_;    
  };
}
