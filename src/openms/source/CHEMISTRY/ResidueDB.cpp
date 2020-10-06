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
    readResiduesFromMap_();
    buildResidueNames_();
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

  void ResidueDB::setResidues_()
  {
    #pragma omp critical (ResidueDB)
    {
      readResiduesFromMap_();
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
    set<String> synonyms = r->getSynonyms();
    for (const String & s : synonyms)
    {
      names.push_back(s);
    }

    if (!r->isModified())
    {
      for (vector<String>::const_iterator it = names.begin(); it != names.end(); ++it)
      {
        residue_names_[*it] = r;
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
      const set<String>& mod_synonyms = mod->getSynonyms();
      for (set<String>::const_iterator it = mod_synonyms.begin(); it != mod_synonyms.end(); ++it)
      {
        mod_names.push_back(*it);
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
    buildResidueNames_();
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

  void ResidueDB::readResiduesFromMap_()
  {
    // create a residue map to avoid the xml file load
    typedef map<String, String> ResidueList;
    ResidueList res;
    
    // fill the ResidueList with the residue information 
    res["Residues:Alanine:Name"] = "Alanine";
    res["Residues:Alanine:ShortName"] =  "Ala";
    res["Residues:Alanine:ThreeLetterCode"] = "ALA";
    res["Residues:Alanine:OneLetterCode"] =  "A";
    res["Residues:Alanine:Formula"] = "C3H7NO2";
    res["Residues:Alanine:pka"] =  "2.35";
    res["Residues:Alanine:pkb"] = "9.87";
    res["Residues:Alanine:GB_SC"] =  "0.00";
    res["Residues:Alanine:GB_BB_L"] = "881.82";
    res["Residues:Alanine:GB_BB_R"] =  "0.00";
    res["Residues:Alanine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Alanine:Synonyms:L-Alanine"] = "L-Alanine";
    res["Residues:Alanine:Synonyms:alanine"] =  "alanine";
    res["Residues:Alanine:Synonyms:Alanin"] = "Alanin";
    res["Residues:Alanine:Synonyms:alanin"] = "alanin";
    res["Residues:Alanine:Synonyms:Ala"] = "Ala";
    
    res["Residues:Cysteine:Name"] =  "Cysteine";
    res["Residues:Cysteine:ShortName"] = "Cys";
    res["Residues:Cysteine:ThreeLetterCode"] = "CYS";
    res["Residues:Cysteine:OneLetterCode"] =  "C";
    res["Residues:Cysteine:Formula"] =  "C3H7NO2S";
    res["Residues:Cysteine:pka"] = "1.92";
    res["Residues:Cysteine:pkb"] = "10.70";
    res["Residues:Cysteine:pkc"] = "8.18";
    res["Residues:Cysteine:GB_SC"] = "0.00";
    res["Residues:Cysteine:GB_BB_L"] = "880.99";
    res["Residues:Cysteine:GB_BB_R"] =  "0.12";
    res["Residues:Cysteine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Cysteine:Synonyms:Cys"] = "Cys";
    res["Residues:Cysteine:Synonyms:Cystine"] = "Cystine";
    
    res["Residues:Aspartate:Name"] =  "Aspartate";
    res["Residues:Aspartate:ShortName"] = "Asp";
    res["Residues:Aspartate:ThreeLetterCode"] = "ASP";
    res["Residues:Aspartate:OneLetterCode"] = "D";
    res["Residues:Aspartate:GB_SC"] = "784.0";
    res["Residues:Aspartate:GB_BB_L"] = "880.02";
    res["Residues:Aspartate:GB_BB_R"] = "-0.63";
    res["Residues:Aspartate:Formula"] = "C4H7NO4";
    res["Residues:Aspartate:pka"] = "1.99";
    res["Residues:Aspartate:pkb"] = "9.90";
    res["Residues:Aspartate:pkc"] = "3.90";
    res["Residues:Aspartate:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Aspartate:Synonyms:Asp"] = "Asp";
    res["Residues:Aspartate:Losses:LossName"] = "water";
    res["Residues:Aspartate:Losses:LossFormula"] = "H2O";
    
    res["Residues:Glutamate:Name" ]= "Glutamate";
    res["Residues:Glutamate:ShortName" ]= "Glu";
    res["Residues:Glutamate:ThreeLetterCode" ]= "GLU";
    res["Residues:Glutamate:OneLetterCode" ]= "E";
    res["Residues:Glutamate:GB_SC" ]= "790.0";
    res["Residues:Glutamate:GB_BB_L" ]= "880.10";
    res["Residues:Glutamate:GB_BB_R" ]= "-0.39";
    res["Residues:Glutamate:Formula" ]= "C5H9NO4";
    res["Residues:Glutamate:pka" ]= "2.10";
    res["Residues:Glutamate:pkb" ]= "9.47";
    res["Residues:Glutamate:pkc" ]= "4.07";
    res["Residues:Glutamate:ResidueSets" ]= "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Glutamate:Synonyms:Glu" ]= "Glu";
    res["Residues:Glutamate:Losses:LossName" ]= "water";
    res["Residues:Glutamate:Losses:LossFormula" ]= "H2O";
    res["Residues:Glutamate:NTermLosses:LossName" ]= "water";
    res["Residues:Glutamate:NTermLosses:LossFormula" ]= "H2O";
    
    res["Residues:Phenylalanine:Name" ]= "Phenylalanine";
    res["Residues:Phenylalanine:ShortName" ]= "Phe";
    res["Residues:Phenylalanine:ThreeLetterCode" ]= "PHE";
    res["Residues:Phenylalanine:OneLetterCode" ]= "F";
    res["Residues:Phenylalanine:Formula" ]= "C9H11NO2";
    res["Residues:Phenylalanine:pka" ]= "2.20";
    res["Residues:Phenylalanine:pkb" ]= "9.31";
    res["Residues:Phenylalanine:GB_SC" ]= "0.00";
    res["Residues:Phenylalanine:GB_BB_L" ]= "881.08";
    res["Residues:Phenylalanine:GB_BB_R" ]= "0.03";
    res["Residues:Phenylalanine:ResidueSets" ]= "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Phenylalanine:Synonyms:Phe" ]= "Phe";
      
    res["Residues:Glycine:Name"] = "Glycine";
    res["Residues:Glycine:ShortName"] = "Gly";
    res["Residues:Glycine:ThreeLetterCode"] = "GLY";
    res["Residues:Glycine:OneLetterCode"] = "G";
    res["Residues:Glycine:Formula"] = "C2H5NO2";
    res["Residues:Glycine:pka"] = "2.35";
    res["Residues:Glycine:pkb"] = "9.78";
    res["Residues:Glycine:GB_SC"] = "0.00";
    res["Residues:Glycine:GB_BB_L"] = "881.17";
    res["Residues:Glycine:GB_BB_R"] = "0.92";
    res["Residues:Glycine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Glycine:Synonyms:Gly"] = "Gly";

    res["Residues:Histidine:Name"] = "Histidine";
    res["Residues:Histidine:ShortName"] = "His";
    res["Residues:Histidine:ThreeLetterCode"] = "HIS";
    res["Residues:Histidine:OneLetterCode"] = "H";
    res["Residues:Histidine:Formula"] = "C6H9N3O2";
    res["Residues:Histidine:pka"] = "1.80";
    res["Residues:Histidine:pkb"] = "9.33";
    res["Residues:Histidine:pkc"] = "6.04";
    res["Residues:Histidine:GB_SC"] = "927.84";
    res["Residues:Histidine:GB_BB_L"] = "881.27";
    res["Residues:Histidine:GB_BB_R"] = "-0.19";
    res["Residues:Histidine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Histidine:Synonyms:His"] = "His";

    res["Residues:Isoleucine:Name"] = "Isoleucine";
    res["Residues:Isoleucine:ShortName"] = "Ile";
    res["Residues:Isoleucine:ThreeLetterCode"] = "ILE";
    res["Residues:Isoleucine:OneLetterCode"] = "I";
    res["Residues:Isoleucine:Formula"] = "C6H13NO2";
    res["Residues:Isoleucine:pka"] = "2.32";
    res["Residues:Isoleucine:pkb"] = "9.76";
    res["Residues:Isoleucine:GB_SC"] = "0.00";
    res["Residues:Isoleucine:GB_BB_L"] = "880.99";
    res["Residues:Isoleucine:GB_BB_R"] = "-1.17";
    res["Residues:Isoleucine:ResidueSets"] = "All,Natural20,Natural19WithoutL,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Isoleucine:Synonyms:Ile"] = "Ile";

    res["Residues:Lysine:Name"] = "Lysine";
    res["Residues:Lysine:ShortName"] = "Lys";
    res["Residues:Lysine:ThreeLetterCode"] = "LYS";
    res["Residues:Lysine:OneLetterCode"] = "K";
    res["Residues:Lysine:Formula"] = "C6H14N2O2";
    res["Residues:Lysine:pka"] = "2.16";
    res["Residues:Lysine:pkb"] = "9.06";
    res["Residues:Lysine:pkc"] = "10.54";
    res["Residues:Lysine:GB_SC"] = "926.74";
    res["Residues:Lysine:GB_BB_L"] = "880.06";
    res["Residues:Lysine:GB_BB_R"] = "-0.71";
    res["Residues:Lysine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Lysine:Synonyms:Lys"] = "Lys";
    res["Residues:Lysine:Losses:LossName"] = "ammonia";
    res["Residues:Lysine:Losses:LossFormula"] = "NH3";

    res["Residues:Leucine:Name"] = "Leucine";
    res["Residues:Leucine:ShortName"] = "Leu";
    res["Residues:Leucine:ThreeLetterCode"] = "LEU";
    res["Residues:Leucine:OneLetterCode"] = "L";
    res["Residues:Leucine:Formula"] = "C6H13NO2";
    res["Residues:Leucine:pka"] = "2.33";
    res["Residues:Leucine:pkb"] = "9.74";
    res["Residues:Leucine:GB_SC"] = "0.00";
    res["Residues:Leucine:GB_BB_L"] = "881.88";
    res["Residues:Leucine:GB_BB_R"] = "-0.09";
    res["Residues:Leucine:ResidueSets"] = "All,Natural20,Natural19WithoutI,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Leucine:Synonyms:Leu"] = "Leu";

    res["Residues:Methionine:Name"] = "Methionine";
    res["Residues:Methionine:ShortName"] = "Met";
    res["Residues:Methionine:ThreeLetterCode"] = "MET";
    res["Residues:Methionine:OneLetterCode"] = "M";
    res["Residues:Methionine:Formula"] = "C5H11NO2S";
    res["Residues:Methionine:pka"] = "2.13";
    res["Residues:Methionine:pkb"] = "9.28";
    res["Residues:Methionine:GB_SC"] = "830.0";
    res["Residues:Methionine:GB_BB_L"] = "881.38";
    res["Residues:Methionine:GB_BB_R"] = "0.30";
    res["Residues:Methionine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Methionine:Synonyms:Met"] = "Met";

    res["Residues:Asparagine:Name"] = "Asparagine";
    res["Residues:Asparagine:ShortName"] = "Asn";
    res["Residues:Asparagine:ThreeLetterCode"] = "ASN";
    res["Residues:Asparagine:OneLetterCode"] = "N";
    res["Residues:Asparagine:GB_SC"] = "864.94";
    res["Residues:Asparagine:GB_BB_L"] = "881.18";
    res["Residues:Asparagine:GB_BB_R"] = "1.56";
    res["Residues:Asparagine:Formula"] = "C4H8N2O3";
    res["Residues:Asparagine:pka"] = "2.14";
    res["Residues:Asparagine:pkb"] = "8.72";
    res["Residues:Asparagine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Asparagine:Synonyms:Asn"] = "Asn";
    res["Residues:Asparagine:Losses:LossName"] = "ammonia";
    res["Residues:Asparagine:Losses:LossFormula"] = "NH3";

    res["Residues:Proline:Name"] = "Proline";
    res["Residues:Proline:ShortName"] = "Pro";
    res["Residues:Proline:ThreeLetterCode"] = "PRO";
    res["Residues:Proline:OneLetterCode"] = "P";
    res["Residues:Proline:Formula"] = "C5H9NO2";
    res["Residues:Proline:pka"] = "1.95";
    res["Residues:Proline:pkb"] = "10.64";
    res["Residues:Proline:GB_SC"] = "0.00";
    res["Residues:Proline:GB_BB_L"] = "881.25";
    res["Residues:Proline:GB_BB_R"] = "11.75";
    res["Residues:Proline:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Proline:Synonyms:Pro"] = "Pro";

    res["Residues:Glutamine:Name"] = "Glutamine";
    res["Residues:Glutamine:ShortName"] = "Gln";
    res["Residues:Glutamine:ThreeLetterCode"] = "GLN";
    res["Residues:Glutamine:OneLetterCode"] = "Q";
    res["Residues:Glutamine:GB_SC"] = "865.25";
    res["Residues:Glutamine:GB_BB_L"] = "881.50";
    res["Residues:Glutamine:GB_BB_R"] = "4.10";
    res["Residues:Glutamine:Formula"] = "C5H10N2O3";
    res["Residues:Glutamine:pka"] = "2.17";
    res["Residues:Glutamine:pkb"] = "9.13";
    res["Residues:Glutamine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Glutamine:Synonyms:Gln"] = "Gln";
    res["Residues:Glutamine:Losses:LossName"] = "ammonia";
    res["Residues:Glutamine:Losses:LossFormula"] = "NH3";
    res["Residues:Glutamine:NTermLosses:LossName"] = "water";
    res["Residues:Glutamine:NTermLosses:LossFormula"] = "H2O";

    res["Residues:Arginine:Name"] = "Arginine";
    res["Residues:Arginine:ShortName"] = "Arg";
    res["Residues:Arginine:ThreeLetterCode"] = "ARG";
    res["Residues:Arginine:OneLetterCode"] = "R";
    res["Residues:Arginine:GB_SC"] = "1000.0";
    res["Residues:Arginine:GB_BB_L"] = "882.98";
    res["Residues:Arginine:GB_BB_R"] = "6.28";
    res["Residues:Arginine:Formula"] = "C6H14N4O2";
    res["Residues:Arginine:pka"] = "1.82";
    res["Residues:Arginine:pkb"] = "8.99";
    res["Residues:Arginine:pkc"] = "12.48";
    res["Residues:Arginine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Arginine:Synonyms:Arg"] = "Arg";
    res["Residues:Arginine:Losses:LossName1"] = "ammonia";
    res["Residues:Arginine:Losses:LossFormula1"] = "NH3";
    res["Residues:Arginine:Losses:LossName2"] = "";
    res["Residues:Arginine:Losses:LossFormula2"] = "NHCNH";
    res["Residues:Arginine:Losses:LossName3"] = "";
    res["Residues:Arginine:Losses:LossFormula3"] = "CONH2";

    res["Residues:Selenocysteine:Name"] = "Selenocysteine";
    res["Residues:Selenocysteine:ShortName"] = "Sec";
    res["Residues:Selenocysteine:ThreeLetterCode"] = "SEC";
    res["Residues:Selenocysteine:OneLetterCode"] = "U";
    res["Residues:Selenocysteine:Formula"] = "C3H7NO2Se";
    res["Residues:Selenocysteine:pkc"] = "5.73";
    res["Residues:Selenocysteine:GB_SC"] = "0.00";
    res["Residues:Selenocysteine:GB_BB_L"] = "880.99";
    res["Residues:Selenocysteine:GB_BB_R"] = "0.12";
    res["Residues:Selenocysteine:ResidueSets"] = "All,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Selenocysteine:Synonyms:Sec"] = "Sec";

    res["Residues:Serine:Name"] = "Serine";
    res["Residues:Serine:ShortName"] = "Ser";
    res["Residues:Serine:ThreeLetterCode"] = "SER";
    res["Residues:Serine:OneLetterCode"] = "S";
    res["Residues:Serine:GB_SC"] = "775.0";
    res["Residues:Serine:GB_BB_L"] = "881.08";
    res["Residues:Serine:GB_BB_R"] = "0.98";
    res["Residues:Serine:Formula"] = "C3H7NO3";
    res["Residues:Serine:pka"] = "2.19";
    res["Residues:Serine:pkb"] = "9.21";
    res["Residues:Serine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Serine:Synonyms:Ser"] = "Ser";
    res["Residues:Serine:Losses:LossName"] = "water";
    res["Residues:Serine:Losses:LossFormula"] = "H2O";

    res["Residues:Threonine:Name"] = "Threonine";
    res["Residues:Threonine:ShortName"] = "Thr";
    res["Residues:Threonine:ThreeLetterCode"] = "THR";
    res["Residues:Threonine:OneLetterCode"] = "T";
    res["Residues:Threonine:Formula"] = "C4H9NO3";
    res["Residues:Threonine:pka"] = "2.09";
    res["Residues:Threonine:pkb"] = "9.10";
    res["Residues:Threonine:GB_SC"] = "780.0";
    res["Residues:Threonine:GB_BB_L"] = "881.14";
    res["Residues:Threonine:GB_BB_R"] = "1.21";
    res["Residues:Threonine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Threonine:Synonyms:Thr"] = "Thr";
    res["Residues:Threonine:Losses:LossName"] = "water";
    res["Residues:Threonine:Losses:LossFormula"] = "H2O";

    res["Residues:Valine:Name"] = "Valine";
    res["Residues:Valine:ShortName"] = "Val";
    res["Residues:Valine:ThreeLetterCode"] = "VAL";
    res["Residues:Valine:OneLetterCode"] = "V";
    res["Residues:Valine:Formula"] = "C5H11NO2";
    res["Residues:Valine:pka"] = "2.39";
    res["Residues:Valine:pkb"] = "9.74";
    res["Residues:Valine:GB_SC"] = "0.0";
    res["Residues:Valine:GB_BB_L"] = "881.17";
    res["Residues:Valine:GB_BB_R"] = "-0.90";
    res["Residues:Valine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Valine:Synonyms:Val"] = "Val";

    res["Residues:Tryptophan:Name"] = "Tryptophan";
    res["Residues:Tryptophan:ShortName"] = "Trp";
    res["Residues:Tryptophan:ThreeLetterCode"] = "TRP";
    res["Residues:Tryptophan:OneLetterCode"] = "W";
    res["Residues:Tryptophan:Formula"] = "C11H12N2O2";
    res["Residues:Tryptophan:pka"] = "2.46";
    res["Residues:Tryptophan:pkb"] = "9.41";
    res["Residues:Tryptophan:GB_SC"] = "909.53";
    res["Residues:Tryptophan:GB_BB_L"] = "881.31";
    res["Residues:Tryptophan:GB_BB_R"] = "0.10";
    res["Residues:Tryptophan:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Tryptophan:Synonyms:Trp"] = "Trp";

    res["Residues:Tyrosine:Name"] = "Tyrosine";
    res["Residues:Tyrosine:ShortName"] = "Tyr";
    res["Residues:Tyrosine:ThreeLetterCode"] = "TYR";
    res["Residues:Tyrosine:OneLetterCode"] = "Y";
    res["Residues:Tyrosine:Formula"] = "C9H11NO3";
    res["Residues:Tyrosine:pka"] = "2.20";
    res["Residues:Tyrosine:pkb"] = "9.21";
    res["Residues:Tyrosine:pkc"] = "10.46";
    res["Residues:Tyrosine:GB_SC"] = "790.0";
    res["Residues:Tyrosine:GB_BB_L"] = "881.20";
    res["Residues:Tyrosine:GB_BB_R"] = "-0.38";
    res["Residues:Tyrosine:ResidueSets"] = "All,Natural20,Natural19WithoutI,Natural19WithoutL,Natural19J,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Tyrosine:Synonyms:Tyr"] = "Tyr";

    res["Residues:Pyrrolysine:Name"] = "Pyrrolysine";
    res["Residues:Pyrrolysine:ShortName"] = "Pyr";
    res["Residues:Pyrrolysine:ThreeLetterCode"] = "PYR";
    res["Residues:Pyrrolysine:OneLetterCode"] = "O";
    res["Residues:Pyrrolysine:Formula"] = "C12H21N3O3";
    res["Residues:Pyrrolysine:pka"] = "0.0";
    res["Residues:Pyrrolysine:pkb"] = "0.0";
    res["Residues:Pyrrolysine:pkc"] = "0.0";
    res["Residues:Pyrrolysine:GB_SC"] = "0.0";
    res["Residues:Pyrrolysine:GB_BB_L"] = "0.0";
    res["Residues:Pyrrolysine:GB_BB_R"] = "0.0";
    res["Residues:Pyrrolysine:ResidueSets"] = "All,AmbiguousWithoutX,Ambiguous,AllNatural";
    res["Residues:Pyrrolysine:Synonyms:Pyr"] = "Pyr";

    res["Residues:Ambiguous_ASN_or_ASP:Name"] = "Asparagine/Aspartate";
    res["Residues:Ambiguous_ASN_or_ASP:ShortName"] = "Asx";
    res["Residues:Ambiguous_ASN_or_ASP:ThreeLetterCode"] = "ASX";
    res["Residues:Ambiguous_ASN_or_ASP:OneLetterCode"] = "B";
    res["Residues:Ambiguous_ASN_or_ASP:Formula"] = "";
    res["Residues:Ambiguous_ASN_or_ASP:GB_SC"] = "0";
    res["Residues:Ambiguous_ASN_or_ASP:GB_BB_L"] = "0";
    res["Residues:Ambiguous_ASN_or_ASP:GB_BB_R"] = "0";
    res["Residues:Ambiguous_ASN_or_ASP:ResidueSets"] = "All,AmbiguousWithoutX,Ambiguous";
    res["Residues:Ambiguous_ASN_or_ASP:Synonyms:Asx"] = "Asx";

    res["Residues:Ambiguous_GLN_or_GLU:Name"] = "Glutamine/Glutamate";
    res["Residues:Ambiguous_GLN_or_GLU:ShortName"] = "Glx";
    res["Residues:Ambiguous_GLN_or_GLU:ThreeLetterCode"] = "GLX";
    res["Residues:Ambiguous_GLN_or_GLU:OneLetterCode"] = "Z";
    res["Residues:Ambiguous_GLN_or_GLU:Formula"] = "";
    res["Residues:Ambiguous_GLN_or_GLU:GB_SC"] = "0";
    res["Residues:Ambiguous_GLN_or_GLU:GB_BB_L"] = "0";
    res["Residues:Ambiguous_GLN_or_GLU:GB_BB_R"] = "0";
    res["Residues:Ambiguous_GLN_or_GLU:ResidueSets"] = "All,AmbiguousWithoutX,Ambiguous";
    res["Residues:Ambiguous_GLN_or_GLU:Synonyms:Glx"] = "Glx";

    res["Residues:Ambiguous_ILE_or_LEU:Name"] = "Isoleucine/Leucine";
    res["Residues:Ambiguous_ILE_or_LEU:ShortName"] = "Xle";
    res["Residues:Ambiguous_ILE_or_LEU:ThreeLetterCode"] = "XLE";
    res["Residues:Ambiguous_ILE_or_LEU:OneLetterCode"] = "J";
    res["Residues:Ambiguous_ILE_or_LEU:Formula"] = "C6H13NO2";
    res["Residues:Ambiguous_ILE_or_LEU:GB_SC"] = "0";
    res["Residues:Ambiguous_ILE_or_LEU:GB_BB_L"] = "880.99";
    res["Residues:Ambiguous_ILE_or_LEU:GB_BB_R"] = "-1.17";
    res["Residues:Ambiguous_ILE_or_LEU:ResidueSets"] = "All,AmbiguousWithoutX,Ambiguous";
    res["Residues:Ambiguous_ILE_or_LEU:Synonyms:Xle"] = "Xle";

    res["Residues:Ambiguous_Unspecified_or_Unknown:Name"] = "Unspecified/Unknown";
    res["Residues:Ambiguous_Unspecified_or_Unknown:ShortName"] = "Xaa";
    res["Residues:Ambiguous_Unspecified_or_Unknown:ThreeLetterCode"] = "XAA";
    res["Residues:Ambiguous_Unspecified_or_Unknown:OneLetterCode"] = "X";
    res["Residues:Ambiguous_Unspecified_or_Unknown:Formula"] = "";
    res["Residues:Ambiguous_Unspecified_or_Unknown:GB_SC"] = "0";
    res["Residues:Ambiguous_Unspecified_or_Unknown:GB_BB_L"] = "0";
    res["Residues:Ambiguous_Unspecified_or_Unknown:GB_BB_R"] = "0";
    res["Residues:Ambiguous_Unspecified_or_Unknown:ResidueSets"] = "All,Ambiguous";
    res["Residues:Ambiguous_Unspecified_or_Unknown:Synonyms:Xaa"] = "Xaa";
    res["Residues:Ambiguous_Unspecified_or_Unknown:Synonyms:Unk"] = "Unk";

    // clear names and lookup
    clearResidues_();
    clearResidueModifications_();

    try
    {
      vector<String> split;
      String first_elem_res = res.begin()->first;
      first_elem_res.split(':', split);
        
      String prefix = split[0] + split[1];
      Residue* res_ptr = nullptr;

      Map<String, String> values;

      for (auto it = res.begin(); it != res.end(); ++it)
      {
        String current_res = it->first;
        current_res.split(':', split);
        
        if (prefix != split[0] + split[1])
        {
          // add residue
          res_ptr = parseResidue_(values);
          values.clear();
          residues_.insert(res_ptr);
          const_residues_.insert(res_ptr);
          prefix = split[0] + split[1];
          residue_by_one_letter_code_[static_cast<unsigned char>(res_ptr->getOneLetterCode()[0])] = res_ptr;
        }
        String value = it->second;
        String key = it->first;
        values[key] = value;
      }
    
      // add last residue
      res_ptr = parseResidue_(values);
      residues_.insert(res_ptr);
      const_residues_.insert(res_ptr);
      residue_by_one_letter_code_[static_cast<unsigned char>(res_ptr->getOneLetterCode()[0])] = res_ptr;
    }
    catch (Exception::BaseException& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, e.what(), "");
    }
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

  Residue* ResidueDB::parseResidue_(Map<String, String>& values)
  {
    vector<EmpiricalFormula> low_mass_ions;
    Residue* res_ptr = new Residue();

    for (Map<String, String>::iterator it = values.begin(); it != values.end(); ++it)
    {
      String key(it->first);
      String value(it->second);

      if (key.hasSuffix(":Name"))
      {
        res_ptr->setName(value);
        continue;
      }
      if (key.hasSuffix(":ShortName"))
      {
        res_ptr->setShortName(value);
        continue;
      }
      if (key.hasSuffix(":ThreeLetterCode"))
      {
        res_ptr->setThreeLetterCode(value);
        continue;
      }
      if (key.hasSuffix(":OneLetterCode"))
      {
        res_ptr->setOneLetterCode(value);
        continue;
      }
      if (key.hasSuffix(":Formula"))
      {
        EmpiricalFormula formula(value);
        res_ptr->setFormula(EmpiricalFormula(value));
        res_ptr->setAverageWeight(formula.getAverageWeight());
        res_ptr->setMonoWeight(formula.getMonoWeight());
        continue;
      }

      if (key.hasSubstring(":Losses:LossName"))
      {
        res_ptr->addLossName(value);
        continue;
      }
      if (key.hasSubstring(":Losses:LossFormula"))
      {
        EmpiricalFormula loss(value);
        res_ptr->addLossFormula(loss);
        continue;
      }

      if (key.hasSubstring("NTermLosses:LossName"))
      {
        res_ptr->addNTermLossName(value);
        continue;
      }

      if (key.hasSubstring("NTermLosses:LossFormula"))
      {
        EmpiricalFormula loss(value);
        res_ptr->addNTermLossFormula(loss);
        continue;
      }

      if (key.hasSubstring("LowMassIons"))
      {
        // no markers defined?
        if (!key.hasSuffix(":"))
        {
          low_mass_ions.push_back(EmpiricalFormula(value));
        }
        continue;
      }
      if (key.hasSubstring("Synonyms"))
      {
        // no synonyms defined?
        if (!key.hasSuffix(":"))
        {
          res_ptr->addSynonym(value);
        }
        continue;
      }
      if (key.hasSubstring("pka"))
      {
        // no pka defined?
        if (!key.hasSuffix(":"))
        {
          res_ptr->setPka(value.toDouble());
        }
        continue;
      }
      if (key.hasSubstring("pkb"))
      {
        // no pkb defined?
        if (!key.hasSuffix(":"))
        {
          res_ptr->setPkb(value.toDouble());
        }
        continue;
      }
      if (key.hasSubstring("pkc"))
      {
        // no pkc defined?
        if (!key.hasSuffix(":"))
        {
          res_ptr->setPkc(value.toDouble());
        }
        continue;
      }
      if (key.hasSubstring("GB_SC"))
      {
        res_ptr->setSideChainBasicity(value.toDouble());
        continue;
      }
      if (key.hasSubstring("GB_BB_L"))
      {
        res_ptr->setBackboneBasicityLeft(value.toDouble());
        continue;
      }
      if (key.hasSubstring("GB_BB_R"))
      {
        res_ptr->setBackboneBasicityRight(value.toDouble());
        continue;
      }
      if (key.hasSubstring("ResidueSets"))
      {
        StringList residue_sets = ListUtils::create<String>(value);
        for (StringList::const_iterator local_it = residue_sets.begin(); local_it != residue_sets.end(); ++local_it)
        {
          res_ptr->addResidueSet(*local_it);
          residue_sets_.insert(*local_it);
        }
        continue;
      }
      cerr << "unknown key: " << key << ", with value: " << value << endl;
    }

    if (!low_mass_ions.empty())
    {
      res_ptr->setLowMassIons(low_mass_ions);
    }

    for (set<String>::const_iterator it = res_ptr->getResidueSets().begin(); it != res_ptr->getResidueSets().end(); ++it)
    {
      residues_by_set_[*it].insert(res_ptr);
    }

    return res_ptr;
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

  void ResidueDB::buildResidueNames_()
  {
    set<Residue*>::iterator it;
    for (it = residues_.begin(); it != residues_.end(); ++it)
    {
      residue_names_[(*it)->getName()] = *it;
      if ((*it)->getThreeLetterCode() != "")
      {
        residue_names_[(*it)->getThreeLetterCode()] = *it;
      }
      if ((*it)->getOneLetterCode() != "")
      {
        residue_names_[(*it)->getOneLetterCode()] = *it;
      }
      if ((*it)->getShortName() != "")
      {
        residue_names_[(*it)->getShortName()] = *it;
      }
      set<String>::iterator sit;
      set<String> syn = (*it)->getSynonyms();
      for (sit = syn.begin(); sit != syn.end(); ++sit)
      {
        if (*sit != "")
        {
          residue_names_[*sit] = *it;
        }
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
