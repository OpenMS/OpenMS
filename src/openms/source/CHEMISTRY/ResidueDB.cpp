// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Jang Jang Jin, Timo Sachsenberg$
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
  ResidueDB::ResidueDB()
  { 
    initResidues_();
  }

  ResidueDB* ResidueDB::getInstance()
  {
    static ResidueDB* db_ = new ResidueDB(); // Meyers' singleton -> thread safe
    return db_;
  }

  ResidueDB::~ResidueDB()
  {
    // free memory
    for (auto& r : const_residues_) { delete r; }
    for (auto& r : const_modified_residues_) { delete r; }
  }

  const Residue* ResidueDB::getResidue(const String& name) const
  {
    if (name.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No residue specified.", "");
    }

    const Residue* r{};
    #pragma omp critical (ResidueDB)
    {   
      auto it = residue_names_.find(name);
      if (it != residue_names_.end()) 
      { 
        r = it->second; 
      }
    }
    if (r == nullptr)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Residue not found: ", name);
    }
    return r;
  }

  const Residue* ResidueDB::getResidue(const unsigned char& one_letter_code) const
  {
    //TODO why does this not throw but the String version does??
    //no lock required here because read only and array is initialized in thread-safe constructor
    return residue_by_one_letter_code_[one_letter_code];
  }

  Size ResidueDB::getNumberOfResidues() const
  {
    Size s;
    #pragma omp critical (ResidueDB)
    {
      s = const_residues_.size();
    } 
    return s;
  }

  Size ResidueDB::getNumberOfModifiedResidues() const
  {
    Size s;
    #pragma omp critical (ResidueDB)
    {
      s = const_modified_residues_.size();
    } 
    return s;
  }

  const set<const Residue*> ResidueDB::getResidues(const String& residue_set) const
  {
    set<const Residue*> s;
    #pragma omp critical (ResidueDB)
    {
      auto it = residues_by_set_.find(residue_set);
      if (it != residues_by_set_.end())
      {
        s = it->second;
      }
    } 

    if (s.empty()) 
    {
      cout << "Residue set cannot be found: '" + residue_set + "'" << endl;
    }
    return s;
  }

  void ResidueDB::initResidues_()
  {
    buildResidues_();    
  }

  void ResidueDB::addResidue_(Residue* r)
  {
    if (!r->isModified())
    { // add (unmodified) residue to residue_names, residues, and const_residues
      const_residues_.insert(r);
      addResidueNames_(r);
    }
    else
    { // add modified residue to const_modified_residues_, and residue_mod_names_
      const_modified_residues_.insert(r);
      addModifiedResidueNames_(r);
    }    
    return;
  }

  bool ResidueDB::hasResidue(const String& res_name) const
  {
    bool found = false;
    #pragma omp critical (ResidueDB)
    {
      found = residue_names_.find(res_name) != residue_names_.end();
    }  
    return found;
  }

  bool ResidueDB::hasResidue(const Residue* residue) const
  {
    bool found = false;
    #pragma omp critical (ResidueDB)
    {
      found = (const_residues_.find(residue) != const_residues_.end() ||
          const_modified_residues_.find(residue) != const_modified_residues_.end());
    } 
    return found;
  }

  void ResidueDB::buildResidues_()
  {
    Residue* alanine = new Residue("Alanine", "Ala", "A", EmpiricalFormula("C3H7NO2"), 2.35, 9.87, -1.00, 0.00, 881.82, 0.00, set<String>{"L-Alanine", "alanine",  "Alanin", "alanin", "ALA"});
    insertResidueAndAssociateWithResidueSet_(alanine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"});
     
    Residue* cysteine = new Residue("Cysteine", "Cys", "C", EmpiricalFormula("C3H7NO2S"), 1.92, 10.70, 8.18, 0.00, 0.12, 880.99,  set<String>{"CYS", "Cystine"});
    insertResidueAndAssociateWithResidueSet_(cysteine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"});
    
    Residue* aspartate = new Residue("Aspartate", "Asp", "D", EmpiricalFormula("C4H7NO4"), 1.99, 9.90, 3.90, 784.0, 880.02, -0.63, set<String>{"ASP"});
    aspartate->addLossName("water");
    aspartate->addLossFormula(EmpiricalFormula("H2O"));
    insertResidueAndAssociateWithResidueSet_(aspartate, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
    
    Residue* glutamate = new Residue ( "Glutamate", "Glu", "E", EmpiricalFormula("C5H9NO4"), 2.10, 9.47, 4.07, 790.0, 880.10, -0.39, set<String>{"GLU"});
    glutamate->addLossName("water");
    glutamate->addLossFormula(EmpiricalFormula("H2O"));
    glutamate->addNTermLossName("water");
    glutamate->addNTermLossFormula(EmpiricalFormula("H2O"));    
    insertResidueAndAssociateWithResidueSet_(glutamate, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
    
    Residue* phenylalanine = new Residue ( "Phenylalanine", "Phe", "F", EmpiricalFormula( "C9H11NO2"), 2.20, 9.31, -1.0, 0.00, 881.08, 0.03, set<String>{"PHE"});
    insertResidueAndAssociateWithResidueSet_(phenylalanine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
    
    Residue* glycine = new Residue ( "Glycine", "Gly", "G", EmpiricalFormula( "C2H5NO2"), 2.35, 9.78, -1.0, 0.00, 881.17, 0.92, set<String>{"GLY"} );
    insertResidueAndAssociateWithResidueSet_(glycine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
    
    Residue* histidin = new Residue ( "Histidine", "His", "H", EmpiricalFormula( "C6H9N3O2"), 1.80, 9.33, 6.04, 927.84, 881.27, -0.19, set<String>{"HIS"});
    insertResidueAndAssociateWithResidueSet_(histidin,  {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"});

    Residue* isoleucine = new Residue ( "Isoleucine", "Ile", "I", EmpiricalFormula( "C6H13NO2"), 2.32, 9.76, -1.0, 0.00, 880.99, -1.17, set<String>{"ILE"});
    insertResidueAndAssociateWithResidueSet_(isoleucine, {"All","Natural20","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"});
      
    Residue* lysine = new Residue ( "Lysine", "Lys", "K", EmpiricalFormula( "C6H14N2O2"), 2.16, 9.06, 10.54, 926.74, 880.06, -0.71, set<String>{ "LYS"});
    lysine->addLossName("ammonia");
    lysine->addLossFormula(EmpiricalFormula("NH3"));
    insertResidueAndAssociateWithResidueSet_(lysine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
	
    Residue* leucine = new Residue ( "Leucine", "Leu", "L", EmpiricalFormula( "C6H13NO2"), 2.33, 9.74, -1.0, 0.00, 881.88, -0.09, set<String>{ "LEU"});
    insertResidueAndAssociateWithResidueSet_(leucine, {"All","Natural20","Natural19WithoutI","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"}  );
	
    Residue* methionine = new Residue ( "Methionine", "Met", "M", EmpiricalFormula( "C5H11NO2S"), 2.13, 9.28, -1.0, 830.0, 881.38, 0.30, set<String>{ "MET"});
    insertResidueAndAssociateWithResidueSet_(methionine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
    
    Residue* asparagine = new Residue ( "Asparagine", "Asn", "N", EmpiricalFormula( "C4H8N2O3"), 2.14, 8.72, -1.0, 864.94, 881.18, 1.56, set<String>{ "ASN"});
    asparagine->addLossName("ammonia");
    asparagine->addLossFormula(EmpiricalFormula("NH3"));
    insertResidueAndAssociateWithResidueSet_(asparagine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"});
	
    Residue* proline = new Residue ( "Proline", "Pro", "P", EmpiricalFormula( "C5H9NO2"), 1.95, 10.64, -1.0, 0.00, 881.25, 11.75, set<String>{ "PRO"});
    insertResidueAndAssociateWithResidueSet_(proline, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
	
    Residue* glutamine = new Residue ( "Glutamine", "Gln", "Q", EmpiricalFormula( "C5H10N2O3"), 2.17, 9.13, -1.0, 865.25, 881.50, 4.10, set<String>{ "GLN"});
    glutamine->addLossName("ammonia");
    glutamine->addLossFormula(EmpiricalFormula("NH3"));
    glutamine->addNTermLossName("water");
    glutamine->addNTermLossFormula(EmpiricalFormula("H2O"));
    insertResidueAndAssociateWithResidueSet_(glutamine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
	
    Residue* arginine = new Residue ( "Arginine", "Arg", "R", EmpiricalFormula( "C6H14N4O2"), 1.82, 8.99, 12.48, 1000.0, 882.98, 6.28, set<String>{ "ARG"});
    arginine->addLossName("ammonia");
    arginine->addLossFormula(EmpiricalFormula("NH3"));
    arginine->addLossName("");
    arginine->addLossFormula(EmpiricalFormula("NHCNH"));
    arginine->addLossName("");
    arginine->addLossFormula(EmpiricalFormula("CONH2"));	  
    insertResidueAndAssociateWithResidueSet_(arginine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
	
    Residue* selenocysteine = new Residue ( "Selenocysteine", "Sec", "U", EmpiricalFormula( "C3H7NO2Se"), 0.00, 0.00, 5.73, 0.00, 880.99, 0.12, set<String>{ "SEC"});
    insertResidueAndAssociateWithResidueSet_(selenocysteine, {"All","AmbiguousWithoutX","Ambiguous","AllNatural"});
    
    Residue* serine = new Residue ( "Serine", "Ser", "S", EmpiricalFormula( "C3H7NO3"), 2.19, 9.21, -1.0, 775.0, 881.08, 0.98, set<String>{ "SER"});
    serine->addLossName("water");
    serine->addLossFormula(EmpiricalFormula("H2O"));
    insertResidueAndAssociateWithResidueSet_(serine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
	
    Residue* threonine = new Residue ( "Threonine", "Thr", "T", EmpiricalFormula( "C4H9NO3"), 2.09, 9.10, -1.0, 780.0, 881.14, 1.21, set<String>{ "THR"});
    threonine->addLossName("water");
    threonine->addLossFormula(EmpiricalFormula("H2O"));
    insertResidueAndAssociateWithResidueSet_(threonine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );
	
    Residue* valine = new Residue ( "Valine", "Val", "V", EmpiricalFormula( "C5H11NO2"), 2.39, 9.74, -1.0, 0.0, 881.17, -0.90, set<String>{ "VAL"});
    insertResidueAndAssociateWithResidueSet_(valine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );

    Residue* tryptophan = new Residue ( "Tryptophan", "Trp", "W", EmpiricalFormula( "C11H12N2O2"), 2.46, 9.41, -1.0, 909.53, 881.31, 0.10, set<String>{ "TRP"});
    insertResidueAndAssociateWithResidueSet_(tryptophan, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );

    Residue* tyrosine = new Residue ( "Tyrosine", "Tyr", "Y", EmpiricalFormula( "C9H11NO3"), 2.20, 9.21, 10.46, 790.0, 881.20, -0.38, set<String>{ "TYR" });
    insertResidueAndAssociateWithResidueSet_(tyrosine, {"All","Natural20","Natural19WithoutI","Natural19WithoutL","Natural19J","AmbiguousWithoutX","Ambiguous","AllNatural"} );

    Residue* pyrrolysine = new Residue ( "Pyrrolysine", "Pyr", "O", EmpiricalFormula( "C12H21N3O3"), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, set<String>{ "PYR"});
    insertResidueAndAssociateWithResidueSet_( pyrrolysine, {"All","AmbiguousWithoutX","Ambiguous","AllNatural"});
	
    Residue* asparagine_aspartate = new Residue ( "Asparagine/Aspartate", "Asx", "B", EmpiricalFormula(""), 0.00, 0.00, -1.0, 0.00, 0.00, 0.00, set<String>{ "ASX" });
    insertResidueAndAssociateWithResidueSet_( asparagine_aspartate , {"All","AmbiguousWithoutX","Ambiguous"});
	
    Residue* glutamine_glutamate = new Residue ( "Glutamine/Glutamate", "Glx", "Z", EmpiricalFormula(""), 0.00, 0.00, -1.0, 0.00, 0.00, 0.00, set<String>{ "GLX"});
    insertResidueAndAssociateWithResidueSet_(glutamine_glutamate, {"All","AmbiguousWithoutX","Ambiguous"} );
	
    Residue* isoleucine_leucine = new Residue ( "Isoleucine/Leucine", "Xle", "J", EmpiricalFormula( "C6H13NO2"), 0.00, 0.00, -1.0, 0.00, 880.99, -1.17, set<String>{ "XLE"});
    insertResidueAndAssociateWithResidueSet_( isoleucine_leucine, {"All","AmbiguousWithoutX","Ambiguous"});

    Residue* unspecified_unknown = new Residue ( "Unspecified/Unknown", "Xaa", "X", EmpiricalFormula(""), 0.00, 0.00, -1.0, 0.00, 0.00, 0.00, set<String>{ "XAA", "Unk"});
    insertResidueAndAssociateWithResidueSet_(unspecified_unknown, {"All","Ambiguous"} );
  }
  
  void ResidueDB::insertResidueAndAssociateWithResidueSet_(Residue* res_ptr, const StringList& residue_sets)
  {    
    for (const String& s : residue_sets)
    {
      res_ptr->addResidueSet(s);
      residue_sets_.insert(s);
    }
        
    for (const String& s : res_ptr->getResidueSets())
    {
      residues_by_set_[s].insert(res_ptr);
    }

    const_residues_.insert(res_ptr);
    residue_by_one_letter_code_[static_cast<unsigned char>(res_ptr->getOneLetterCode()[0])] = res_ptr;

    addResidueNames_(res_ptr);
   }

  const set<String> ResidueDB::getResidueSets() const
  {
    set<String> rs;
    #pragma omp critical (ResidueDB)
    {
      rs = residue_sets_; 
    }
    return rs;
  }

  void ResidueDB::addModifiedResidueNames_(const Residue* r)
  {
    // get all modification names
    vector<String> mod_names;
    const ResidueModification* mod = r->getModification();

    if (!mod->getId().empty()) mod_names.push_back(mod->getId()); // for user-defined mods
    mod_names.push_back(mod->getFullName());
    mod_names.push_back(mod->getFullId());

    for (const String& s : mod->getSynonyms())
    {
      mod_names.push_back(s);
    }

    vector<String> names;
    // add name to lookup
    if (!r->getName().empty()) 
    {
      names.push_back(r->getName());
    }
    // add all synonyms to lookup
    for (const String & s : r->getSynonyms())
    {
      names.push_back(s);
    }

    for (const String& n : names)
    {
      if (n.empty()) continue;
      for (const String& m : mod_names)
      {
        if (m.empty()) continue;
        residue_mod_names_[n][m] = r;
      }
    }
  }

  void ResidueDB::addResidueNames_(const Residue* r)
  {
    // add name to residue_names_
    residue_names_[r->getName()] = r;

    // add tree letter code to residue_names_
    if (!r->getThreeLetterCode().empty())
    {
      residue_names_[r->getThreeLetterCode()] = r;
    }

    // add one letter code to residue_names_
    if (!r->getOneLetterCode().empty())
    {
      residue_names_[r->getOneLetterCode()] = r;
    }

    // add all synonyms to residue_names_
    for (const String& s : r->getSynonyms())
    {
      if (!s.empty())
      {
        residue_names_[s] = r;
      }
    }
  }

  const Residue* ResidueDB::getModifiedResidue(const String& modification)
  {
    // throws if modification is not part of ModificationsDB
    const ResidueModification* mod = ModificationsDB::getInstance()->getModification(modification, "", ResidueModification::ANYWHERE);
    auto r = getResidue(mod->getOrigin());
    return getModifiedResidue(r, mod->getFullId());
  }

  const Residue* ResidueDB::getModifiedResidue(const Residue* residue, const String& modification)
  {
    OPENMS_PRECONDITION(!modification.empty(), "Modification cannot be empty")
    // search if the mod already exists
    const String & res_name = residue->getName();
    Residue* res{};
    bool residue_found(true), mod_found(true);
    #pragma omp critical (ResidueDB)
    {
      // Perform a single lookup of the residue name in our database, we assume
      // that if it is present in residue_mod_names_ then we have seen it
      // before and can directly grab it. If its not present, we may have as
      // unmodified residue in residue_names_ but need to create a new entry as
      // modified residue. If the residue itself is unknown, we will throw (see
      // below).
      const auto& rm_entry = residue_mod_names_.find(res_name);
      if (rm_entry == residue_mod_names_.end())
      {
        if (residue_names_.find(res_name) == residue_names_.end())
        {
          residue_found = false;
        }
      }

      if (residue_found)
      {
        const ResidueModification* mod{};
                  
        try
        {
          static const ModificationsDB* mdb = ModificationsDB::getInstance();
          if (modification.hasSubstring("-term "))
          {
            // handle terminal modifications of format: "MOD_NAME (Protein {N|C}-term RESIDUE_NAME)"
            if (modification.hasSubstring("Protein N-term"))
            {
              mod = mdb->getModification(modification, residue->getOneLetterCode(), ResidueModification::PROTEIN_N_TERM); 
            } 
            else if (modification.hasSubstring("Protein C-term"))
            {
              mod = mdb->getModification(modification, residue->getOneLetterCode(), ResidueModification::PROTEIN_C_TERM); 
            }
            // handle terminal modifications of format: "MOD_NAME ({N|C}-term RESIDUE_NAME)"
            else if (modification.hasSubstring("N-term"))
            {
              mod = mdb->getModification(modification, residue->getOneLetterCode(), ResidueModification::N_TERM); 
            } 
            else if (modification.hasSubstring("C-term"))
            {
              mod = mdb->getModification(modification, residue->getOneLetterCode(), ResidueModification::C_TERM); 
            }
          }
          else
          {
            mod = mdb->getModification(modification, residue->getOneLetterCode(), ResidueModification::ANYWHERE);
          }  
        }
        catch (...)
        {
          mod_found = false;
        }

        // check if modification in ResidueDB
        if (mod_found)
        {
          // check if modified residue is already present in ResidueDB
          bool found = false;
          if (rm_entry != residue_mod_names_.end())
          {
            const String& id = mod->getId().empty() ? mod->getFullId() : mod->getId();
            const auto& inner = rm_entry->second.find(id);
            if (inner != rm_entry->second.end())
            {
              res = const_cast<Residue*>(inner->second);
              found = true;
            }
          }
          if (!found)
          {
            // create and register this modified residue
            res = new Residue(*residue_names_.at(res_name));
            res->setModification(mod);
            addResidue_(res);            
          }
        }
      }
    }

    // throwing (uncaught) exceptions needs to happen outside of critical section (see OpenMP reference manual)
    if (!residue_found)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Residue not found: ", res_name);
    }
    else if (!mod_found)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Modification not found: ", modification);
    }
    
    return res;
  }

  const Residue* ResidueDB::getModifiedResidue(const Residue* residue, const ResidueModification* mod)
  {
    OPENMS_PRECONDITION(mod != nullptr, "Mod cannot be nullptr")
    OPENMS_PRECONDITION(mod->getTermSpecificity() == ResidueModification::ANYWHERE, "Mod's term specificity needs to be ANYWHERE to attach it to Residues");
    OPENMS_PRECONDITION(mod->getOrigin() == residue->getOneLetterCode()[0], "Mod's AA origin needs to match residues one-letter-code");
    // search if the mod already exists
    const String & res_name = residue->getName();
    Residue* res{};
    bool residue_found(true);
    #pragma omp critical (ResidueDB)
    {
      // Perform a single lookup of the residue name in our database, we assume
      // that if it is present in residue_mod_names_ then we have seen it
      // before and can directly grab it. If its not present, we may have as
      // unmodified residue in residue_names_ but need to create a new entry as
      // modified residue. If the residue itself is unknown, we will throw (see
      // below).
      const auto& rm_entry = residue_mod_names_.find(res_name);
      if (rm_entry == residue_mod_names_.end())
      {
        if (residue_names_.find(res_name) == residue_names_.end())
        {
          residue_found = false;
        }
      }

      if (residue_found)
      {

        // check if modification in ResidueDB
        if (mod != nullptr)
        {
          // check if modified residue is already present in ResidueDB
          bool found = false;
          if (rm_entry != residue_mod_names_.end())
          {
            const String& id = mod->getId().empty() ? mod->getFullId() : mod->getId();
            const auto& inner = rm_entry->second.find(id);
            if (inner != rm_entry->second.end())
            {
              res = const_cast<Residue*>(inner->second);
              found = true;
            }
          }
          if (!found)
          {
            // create and register this modified residue
            res = new Residue(*residue_names_.at(res_name));
            res->setModification(mod);
            addResidue_(res);            
          }
        }
      }
    }

    // throwing (uncaught) exceptions needs to happen outside of critical section (see OpenMP reference manual)
    if (!residue_found)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Residue not found: ", res_name);
    }
    
    return res;
  }
}
