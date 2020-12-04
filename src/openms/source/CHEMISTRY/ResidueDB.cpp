// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
//

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/SYSTEM/File.h>

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
    static ResidueDB* db_ = new ResidueDB;
    return db_;
  }

  ResidueDB::~ResidueDB()
  {
    clear_();
  }

  const Residue* ResidueDB::getResidue(const String& name) const
  {
    if (name.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No residue specified.", "");
    }

    Residue* r(nullptr);
    #pragma omp critical (ResidueDB)
    {   
      auto it = residue_names_.find(name);
      if (it != residue_names_.end()) r = it->second;
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
      s = residues_.size();
    } 
    return s;
  }

  Size ResidueDB::getNumberOfModifiedResidues() const
  {
    Size s;
    #pragma omp critical (ResidueDB)
    {
      s = modified_residues_.size();
    } 
    return s;
  }

  const set<const Residue*> ResidueDB::getResidues(const String& residue_set) const
  {
    set<const Residue*> s;
    #pragma omp critical (ResidueDB)
    {
      if (residues_by_set_.has(residue_set))
      {
        s = residues_by_set_[residue_set];
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
    #pragma omp critical (ResidueDB)
    {
      buildResidues_();
      buildResidueNames_();
    }     
  }

  void ResidueDB::addResidue_(Residue* r)
  {
    vector<String> names;
    if (r->getName() != "")
    {
      names.push_back(r->getName());
    }
    if (r->getShortName() != "")
    {
      names.push_back(r->getShortName());
    }
    const set<String>& synonyms = r->getSynonyms();
    for (const String & s : synonyms)
    {
      names.push_back(s);
    }

    if (!r->isModified())
    {
      for (const String& name : names)
      {
        residue_names_[name] = r;
      }
      residues_.insert(r);
      const_residues_.insert(r);
    }
    else
    {
      modified_residues_.insert(r);
      const_modified_residues_.insert(r);

      // get all modification names
      vector<String> mod_names;
      const ResidueModification* mod = r->getModification();

      mod_names.push_back(mod->getId());
      mod_names.push_back(mod->getFullName());
      mod_names.push_back(mod->getFullId());

      for (const String& s : mod->getSynonyms())
      {
        mod_names.push_back(s);
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
    buildResidueName_(r);
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
    // clear names and lookup
    clear_();

    Residue* alanine = new  Residue("Alanine", "Ala", "ALA", "A", EmpiricalFormula("C3H7NO2"), EmpiricalFormula("C3H7NO2").getAverageWeight(), EmpiricalFormula("C3H7NO2").getMonoWeight(), 2.35, 9.87, -1.00, 0.00, 881.82, 0.00, set<String>{"L-Alanine", "alanine",  "Alanin", "alanin", "Ala"});
    insertResidues_(alanine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural");
     
    Residue* cysteine = new Residue("Cysteine", "Cys", "CYS", "C", EmpiricalFormula("C3H7NO2S"), EmpiricalFormula("C3H7NO2S").getAverageWeight(), EmpiricalFormula("C3H7NO2S").getMonoWeight(), 1.92, 10.70, 8.18, 0.00, 0.12, 880.99,  set<String>{"Cys", "Cystine"});
    insertResidues_(cysteine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural");
    
    Residue* aspartate = new Residue("Aspartate", "Asp", "ASP", "D", EmpiricalFormula("C4H7NO4"), EmpiricalFormula("C4H7NO4").getAverageWeight(), EmpiricalFormula("C4H7NO4").getMonoWeight(), 1.99, 9.90, 3.90, 784.0, 880.02, -0.63, set<String>{"Asp"});
    aspartate->addLossName("water");
    aspartate->addLossFormula(EmpiricalFormula("H2O"));
    insertResidues_(aspartate, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
    
    Residue* glutamate = new Residue ( "Glutamate", "Glu", "GLU", "E", EmpiricalFormula( "C5H9NO4"), EmpiricalFormula( "C5H9NO4").getAverageWeight(), EmpiricalFormula( "C5H9NO4").getMonoWeight(), 2.10, 9.47, 4.07, 790.0, 880.10, -0.39, set<String>{"Glu"});
    glutamate->addLossName("water");
    glutamate->addLossFormula(EmpiricalFormula("H2O"));
    glutamate->addNTermLossName("water");
    glutamate->addNTermLossFormula(EmpiricalFormula("H2O"));
    insertResidues_(glutamate, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
    
    Residue* phenylalanine = new Residue ( "Phenylalanine", "Phe", "PHE", "F", EmpiricalFormula( "C9H11NO2"), EmpiricalFormula( "C9H11NO2").getAverageWeight(), EmpiricalFormula( "C9H11NO2").getMonoWeight(), 2.20, 9.31, -1.0, 0.00, 881.08, 0.03, set<String>{"Phe"});
  	insertResidues_(phenylalanine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
    
    Residue* glycine = new Residue ( "Glycine", "Gly", "GLY", "G", EmpiricalFormula( "C2H5NO2"), EmpiricalFormula( "C2H5NO2").getAverageWeight(), EmpiricalFormula( "C2H5NO2").getMonoWeight(), 2.35, 9.78, -1.0, 0.00, 881.17, 0.92, set<String>{"Gly"} );
	  insertResidues_(glycine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
    
    Residue* histidin = new Residue ( "Histidine", "His", "HIS", "H", EmpiricalFormula( "C6H9N3O2"), EmpiricalFormula( "C6H9N3O2").getAverageWeight(), EmpiricalFormula( "C6H9N3O2").getMonoWeight(), 1.80, 9.33, 6.04, 927.84, 881.27, -0.19, set<String>{"His"});
	  insertResidues_(histidin,  "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural");

    Residue* isoleucine = new Residue ( "Isoleucine", "Ile", "ILE", "I", EmpiricalFormula( "C6H13NO2"), EmpiricalFormula( "C6H13NO2").getAverageWeight(), EmpiricalFormula( "C6H13NO2").getMonoWeight(), 2.32, 9.76, -1.0, 0.00, 880.99, -1.17, set<String>{"Ile"});
    insertResidues_(isoleucine, "All,Natural20,Natural19WithoutL,AmbiguousWithoutX,Ambiguous,AllNatural" );
      
    Residue* lysine = new Residue ( "Lysine", "Lys", "LYS", "K", EmpiricalFormula( "C6H14N2O2"), EmpiricalFormula( "C6H14N2O2").getAverageWeight(), EmpiricalFormula( "C6H14N2O2").getMonoWeight(), 2.16, 9.06, 10.54, 926.74, 880.06, -0.71, set<String>{ "Lys"});
    lysine->addLossName("ammonia");
    lysine->addLossFormula(EmpiricalFormula("NH3"));
	  insertResidues_(lysine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
	
	  Residue* leucine = new Residue ( "Leucine", "Leu", "LEU", "L", EmpiricalFormula( "C6H13NO2"), EmpiricalFormula( "C6H13NO2").getAverageWeight(), EmpiricalFormula( "C6H13NO2").getMonoWeight(), 2.33, 9.74, -1.0, 0.00, 881.88, -0.09, set<String>{ "Leu"});
	  insertResidues_(leucine, "All,Natural20,Natural19WithoutI,AmbiguousWithoutX,Ambiguous,AllNatural"  );
	
    Residue* methionine = new Residue ( "Methionine", "Met", "MET", "M", EmpiricalFormula( "C5H11NO2S"), EmpiricalFormula( "C5H11NO2S").getAverageWeight(), EmpiricalFormula( "C5H11NO2S").getMonoWeight(), 2.13, 9.28, -1.0, 830.0, 881.38, 0.30, set<String>{ "Met"});
    insertResidues_(methionine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
    
    Residue* asparagine = new Residue ( "Asparagine", "Asn", "ASN", "N", EmpiricalFormula( "C4H8N2O3"), EmpiricalFormula( "C4H8N2O3").getAverageWeight(), EmpiricalFormula( "C4H8N2O3").getMonoWeight(), 2.14, 8.72, -1.0, 864.94, 881.18, 1.56, set<String>{ "Asn"});
    asparagine->addLossName("ammonia");
    asparagine->addLossFormula(EmpiricalFormula("NH3"));
  	insertResidues_(asparagine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural");
	
    Residue* proline = new Residue ( "Proline", "Pro", "PRO", "P", EmpiricalFormula( "C5H9NO2"), EmpiricalFormula( "C5H9NO2").getAverageWeight(), EmpiricalFormula( "C5H9NO2").getMonoWeight(), 1.95, 10.64, -1.0, 0.00, 881.25, 11.75, set<String>{ "Pro"});
    insertResidues_(proline, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
	
    Residue* glutamine = new Residue ( "Glutamine", "Gln", "GLN", "Q", EmpiricalFormula( "C5H10N2O3"), EmpiricalFormula( "C5H10N2O3").getAverageWeight(), EmpiricalFormula( "C5H10N2O3").getMonoWeight(), 2.17, 9.13, -1.0, 865.25, 881.50, 4.10, set<String>{ "Gln"});
    glutamine->addLossName("ammonia");
    glutamine->addLossFormula(EmpiricalFormula("NH3"));
    glutamine->addNTermLossName("water");
    glutamine->addNTermLossFormula(EmpiricalFormula("H2O"));
  	insertResidues_(glutamine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
	
    Residue* arginine = new Residue ( "Arginine", "Arg", "ARG", "R", EmpiricalFormula( "C6H14N4O2"), EmpiricalFormula( "C6H14N4O2").getAverageWeight(), EmpiricalFormula( "C6H14N4O2").getMonoWeight(), 1.82, 8.99, 12.48, 1000.0, 882.98, 6.28, set<String>{ "Arg"});
    arginine->addLossName("ammonia");
    arginine->addLossFormula(EmpiricalFormula("NH3"));
    arginine->addLossName("");
    arginine->addLossFormula(EmpiricalFormula("NHCNH"));
    arginine->addLossName("");
    arginine->addLossFormula(EmpiricalFormula("CONH2"));
	  insertResidues_(arginine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
	
	  Residue* selenocysteine = new Residue ( "Selenocysteine", "Sec", "SEC", "U", EmpiricalFormula( "C3H7NO2Se"), EmpiricalFormula( "C3H7NO2Se").getAverageWeight(), EmpiricalFormula( "C3H7NO2Se").getMonoWeight(), 0.00, 0.00, 5.73, 0.00, 880.99, 0.12, set<String>{ "Sec"});
    insertResidues_(selenocysteine, "All,AmbiguousWithoutX,Ambiguous,AllNatural" );
    
    Residue* serine = new Residue ( "Serine", "Ser", "SER", "S", EmpiricalFormula( "C3H7NO3"), EmpiricalFormula( "C3H7NO3").getAverageWeight(), EmpiricalFormula( "C3H7NO3").getMonoWeight(), 2.19, 9.21, -1.0, 775.0, 881.08, 0.98, set<String>{ "Ser"});
    serine->addLossName("water");
    serine->addLossFormula(EmpiricalFormula("H2O"));
  	insertResidues_(serine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
	
    Residue* threonine = new Residue ( "Threonine", "Thr", "THR", "T", EmpiricalFormula( "C4H9NO3"), EmpiricalFormula( "C4H9NO3").getAverageWeight(), EmpiricalFormula( "C4H9NO3").getMonoWeight(), 2.09, 9.10, -1.0, 780.0, 881.14, 1.21, set<String>{ "Thr"});
  	threonine->addLossName("water");
    threonine->addLossFormula(EmpiricalFormula("H2O"));
	  insertResidues_(threonine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );
	
    Residue* valine = new Residue ( "Valine", "Val", "VAL", "V", EmpiricalFormula( "C5H11NO2"), EmpiricalFormula( "C5H11NO2").getAverageWeight(), EmpiricalFormula( "C5H11NO2").getMonoWeight(), 2.39, 9.74, -1.0, 0.0, 881.17, -0.90, set<String>{ "Val"});
    insertResidues_(valine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );

    Residue* tryptophan = new Residue ( "Tryptophan", "Trp", "TRP", "W", EmpiricalFormula( "C11H12N2O2"), EmpiricalFormula( "C11H12N2O2").getAverageWeight(), EmpiricalFormula( "C11H12N2O2").getMonoWeight(), 2.46, 9.41, -1.0, 909.53, 881.31, 0.10, set<String>{ "Trp"});
    insertResidues_(tryptophan, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );

    Residue* tyrosine = new Residue ( "Tyrosine", "Tyr", "TYR", "Y", EmpiricalFormula( "C9H11NO3"), EmpiricalFormula( "C9H11NO3").getAverageWeight(), EmpiricalFormula( "C9H11NO3").getMonoWeight(), 2.20, 9.21, 10.46, 790.0, 881.20, -0.38, set<String>{ "Tyr" });
    insertResidues_(tyrosine, "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural" );

    Residue* pyrrolysine = new Residue ( "Pyrrolysine", "Pyr", "PYR", "O", EmpiricalFormula( "C12H21N3O3"), EmpiricalFormula( "C12H21N3O3").getAverageWeight(), EmpiricalFormula( "C12H21N3O3").getMonoWeight(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, set<String>{ "Pyr"});
    insertResidues_( pyrrolysine, "All,AmbiguousWithoutX,Ambiguous,AllNatural");
	
    Residue* asparagine_aspartate = new Residue ( "Asparagine/Aspartate", "Asx", "ASX", "B", EmpiricalFormula( ""), EmpiricalFormula( "").getAverageWeight(), EmpiricalFormula( "").getMonoWeight(), 0.00, 0.00, -1.0, 0.00, 0.00, 0.00, set<String>{ "Asx" });
    insertResidues_( asparagine_aspartate , "All,AmbiguousWithoutX,Ambiguous");
	
    Residue* glutamine_glutamate = new Residue ( "Glutamine/Glutamate", "Glx", "GLX", "Z", EmpiricalFormula( ""), EmpiricalFormula( "").getAverageWeight(), EmpiricalFormula( "").getMonoWeight(), 0.00, 0.00, -1.0, 0.00, 0.00, 0.00, set<String>{ "Glx"});
    insertResidues_(glutamine_glutamate, "All,AmbiguousWithoutX,Ambiguous" );
	
    Residue* isoleucine_leucine = new Residue ( "Isoleucine/Leucine", "Xle", "XLE", "J", EmpiricalFormula( "C6H13NO2"), EmpiricalFormula( "C6H13NO2").getAverageWeight(), EmpiricalFormula( "C6H13NO2").getMonoWeight(), 0.00, 0.00, -1.0, 0.00, 880.99, -1.17, set<String>{ "Xle"});
    insertResidues_( isoleucine_leucine, "All,AmbiguousWithoutX,Ambiguous");

    Residue* unspecified_unknown = new Residue ( "Unspecified/Unknown", "Xaa", "XAA", "X", EmpiricalFormula( ""), EmpiricalFormula( "").getAverageWeight(), EmpiricalFormula( "").getMonoWeight(), 0.00, 0.00, -1.0, 0.00, 0.00, 0.00, set<String>{ "Xaa", "Unk"});
    insertResidues_(unspecified_unknown, "All,Ambiguous" );
  }
  
  void ResidueDB::clear_()
  {
    clearResidues_();
    clearResidueModifications_();
  }

  void ResidueDB::clearResidues_()
  {
    // initialize lookup table to null pointer
    for (Size i = 0; i != sizeof(residue_by_one_letter_code_)/sizeof(residue_by_one_letter_code_[0]); ++i)
    {
      residue_by_one_letter_code_[i] = nullptr;
    }

    for (auto& r : residues_) { delete r; }
    residues_.clear();
    residue_names_.clear();
    const_residues_.clear();
    residues_by_set_.clear();
    residue_sets_.clear();
  }

  void ResidueDB::clearResidueModifications_()
  {
    for (auto& r : modified_residues_) { delete r; }
    modified_residues_.clear();
    residue_mod_names_.clear();
    const_modified_residues_.clear();
  }
  
  void ResidueDB::insertResidues_(Residue* res_ptr, const String& residue_set_names)
  {
     StringList residue_sets = ListUtils::create<String>(residue_set_names);
    
    for (const String& s : residue_sets)
        {
          res_ptr->addResidueSet(s);
          residue_sets_.insert(s);
        }
        
   for (const String& s : res_ptr->getResidueSets())
    {
      residues_by_set_[s].insert(res_ptr);
    }

     residues_.insert(res_ptr);
     const_residues_.insert(res_ptr);
     residue_by_one_letter_code_[static_cast<unsigned char>(res_ptr->getOneLetterCode()[0])] = res_ptr;
   
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

  void ResidueDB::buildResidueName_(Residue* r)
  {
    // add name to residue_names_
    residue_names_[r->getName()] = r;

    // add tree letter code to residue_names_
    if (r->getThreeLetterCode() != "")
    {
      residue_names_[r->getThreeLetterCode()] = r;
    }

    // add one letter code to residue_names_
    if (r->getOneLetterCode() != "")
    {
      residue_names_[r->getOneLetterCode()] = r;
    }

    // add short name to residue_names_
    if (r->getShortName() != "")
    {
      residue_names_[r->getShortName()] = r;
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

  void ResidueDB::buildResidueNames_()
  {
    for (Residue* r : residues_)
    {
      buildResidueName_(r)
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
    Residue* res(nullptr);
    bool residue_found(true), mod_found(true);
    #pragma omp critical (ResidueDB)
    {
      // Perform a single lookup of the residue name in our database, we assume
      // that if it is present in residue_mod_names_ then we have seen it
      // before and can directly grab it. If its not present, we may have as
      // unmodified residue in residue_names_ but need to create a new entry as
      // modified residue. If the residue itself is unknow, we will throw (see
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
        const ResidueModification* mod;
        try
        {
          // terminal modifications don't apply to residues (side chain), so only consider internal ones
          static const ModificationsDB* mdb = ModificationsDB::getInstance();
          mod = mdb->getModification(modification, residue->getOneLetterCode(), ResidueModification::ANYWHERE);
        }
        catch (...)
        {
          mod_found = false;
        }

        // check if modification in ResidueDB
        if (mod_found)
        {
          const String& id = mod->getId().empty() ? mod->getFullId() : mod->getId();

          // check if modified residue is already present in ResidueDB
          bool found = false;
          if (rm_entry != residue_mod_names_.end())
          {
            const auto& inner = rm_entry->second.find(id);
            if (inner != rm_entry->second.end())
            {
              res = inner->second;
              found = true;
            }
          }
          if (!found)
          {
            // create and register this modified residue
            res = new Residue(*residue_names_[res_name]);
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

}
