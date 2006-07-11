// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------
// $Maintainer:$
// --------------------------------------------------------------------------



#include "AminoAcid.h"


//STL includes


//My includes
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/FORMAT/Param.h>


using namespace std;
using namespace OpenMS;


void AminoAcid::getID_(string type) throw()
{
  //execute MySQL-query
  sql_adapter_->executeQuery("SELECT aminoacid_ID FROM " + (string)AMINO_TABLE + " WHERE one_letter_code   = \"" + type +
                             "\" OR three_letter_code    = \"" + type +
                             "\" OR    aminoacid_name    = \"" + type + "\"");

  string ID_string;
  if(sql_adapter_->ifGetUnaryResult(ID_string))
    {
      ID = atoi(ID_string.c_str());
    }
  else
    {
      throw UnknownAminoAcid(__FILE__, __LINE__, __PRETTY_FUNCTION__,type);
    }
}


void AminoAcid::initialize_(string type)
{
#ifdef ANNOTATE_XML
  //now we do not have a database
  //we need the name of the aminoacid to identify it
  name_ = code_names_[type];

  Param param;
  param.load( XML_FILE );
  middle_formula_ = (string)param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":middle_formula");
  n_term_formula_ = (string)param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":n_term_formula");
  c_term_formula_ = (string)param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":c_term_formula");

  middle_mono_mass_ = param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":middle_mono_mass");
  n_term_mono_mass_ = param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":n_term_mono_mass");
  c_term_mono_mass_ = param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":c_term_mono_mass");

  middle_average_mass_ = param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":middle_average_mass");
  n_term_average_mass_ = param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":n_term_average_mass");
  c_term_average_mass_ = param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":c_term_average_mass");

  one_letter_code_ = (string)param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":one_letter_code");
  three_letter_code_ = (string)param.getValue("Preferences:SpecAnnotate:Aminoacid:" + name_ + ":three_letter_code");

#endif

#ifndef ANNOTATE_XML
  //now we have a database
  //connect to database
  sql_adapter_ = new MySQLAdapter();
  sql_adapter_->connect(db_username_.c_str(), db_password_.c_str(), db_host_.c_str());
  sql_adapter_->selectDB(DATABASE);

  //get unique ID from for this instance of AminoAcid from database
  getID_(type);

  //get molecular formulae:
  sql_adapter_->executeQuery("SELECT middle_formula FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  middle_formula_ = sql_adapter_->getUnaryResult();
  sql_adapter_->executeQuery("SELECT n_term_formula FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  n_term_formula_ = sql_adapter_->getUnaryResult();
  sql_adapter_->executeQuery("SELECT c_term_formula FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  c_term_formula_ = sql_adapter_->getUnaryResult();
  sql_adapter_->executeQuery("SELECT single_formula FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  single_formula_ = sql_adapter_->getUnaryResult();

  //get monoisotopic masses:
  sql_adapter_->executeQuery("SELECT middle_mono_mass FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  middle_mono_mass_ = atof((sql_adapter_->getUnaryResult()).c_str());
  sql_adapter_->executeQuery("SELECT n_term_mono_mass FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  n_term_mono_mass_ = atof((sql_adapter_->getUnaryResult()).c_str());
  sql_adapter_->executeQuery("SELECT c_term_mono_mass FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  c_term_mono_mass_ = atof((sql_adapter_->getUnaryResult()).c_str());
  sql_adapter_->executeQuery("SELECT single_mono_mass FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  single_mono_mass_ = atof((sql_adapter_->getUnaryResult()).c_str());

  //get average masses:
  sql_adapter_->executeQuery("SELECT middle_average_mass FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  middle_average_mass_ = atof((sql_adapter_->getUnaryResult()).c_str());
  sql_adapter_->executeQuery("SELECT n_term_average_mass FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  n_term_average_mass_ = atof((sql_adapter_->getUnaryResult()).c_str());
  sql_adapter_->executeQuery("SELECT c_term_average_mass FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  c_term_average_mass_ = atof((sql_adapter_->getUnaryResult()).c_str());
  sql_adapter_->executeQuery("SELECT single_average_mass FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  single_average_mass_ = atof((sql_adapter_->getUnaryResult()).c_str());

  //get aminoacid name
  sql_adapter_->executeQuery("SELECT aminoacid_name FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  name_ = sql_adapter_->getUnaryResult();

  //get one-letter-code
  sql_adapter_->executeQuery("SELECT one_letter_code FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  one_letter_code_ = sql_adapter_->getUnaryResult();

  //get three-letter-code
  sql_adapter_->executeQuery("SELECT three_letter_code FROM " + (string)AMINO_TABLE + " WHERE aminoacid_id = \"" + String(ID) + "\"");
  three_letter_code_ = sql_adapter_->getUnaryResult();
#endif
}


AminoAcid::AminoAcid()
{
  cerr << "If you use the non-detailed constructor of class AminoAcid, you should know what you are doing!" << endl;
}


AminoAcid::AminoAcid(string type, string db_username, string db_password, string db_host)
{
#ifdef ANNOTATE_XML
  //initialize the mapping from one-letter-code and three letter code to aminoacid names (only needed for non-database mode)
  code_names_["A"]="Alanine";
  code_names_["ALA"]="Alanine";
  code_names_["Alanine"]="Alanine";
  code_names_["R"]="Arginine";
  code_names_["ARG"]="Arginine";
  code_names_["Arginine"]="Arginine";
  code_names_["N"]="Asparagine";
  code_names_["ASN"]="Asparagine";
  code_names_["Asparagine"]="Asparagine";
  code_names_["D"]="AsparticAcid";
  code_names_["ASP"]="AsparticAcid";
  code_names_["AsparticAcid"]="AsparticAcid";
  code_names_["C"]="Cysteine";
  code_names_["CYS"]="Cysteine";
  code_names_["Cysteine"]="Cysteine";
  code_names_["E"]="GlutamicAcid";
  code_names_["GLU"]="GlutamicAcid";
  code_names_["GlutamicAcid"]="GlutamicAcid";
  code_names_["Q"]="Glutamine";
  code_names_["GLN"]="Glutamine";
  code_names_["Glutamine"]="Glutamine";
  code_names_["G"]="Glycine";
  code_names_["GLY"]="Glycine";
  code_names_["Glycine"]="Glycine";
  code_names_["H"]="Histidine";
  code_names_["HIS"]="Histidine";
  code_names_["Histidine"]="Histidine";
  code_names_["I"]="Isoleucine";
  code_names_["ILE"]="Isoleucine";
  code_names_["Isoleucine"]="Isoleucine";
  code_names_["L"]="Leucine";
  code_names_["LEU"]="Leucine";
  code_names_["Leucine"]="Leucine";
  code_names_["K"]="Lysine";
  code_names_["LYS"]="Lysine";
  code_names_["Lysine"]="Lysine";
  code_names_["M"]="Methionine";
  code_names_["MET"]="Methionine";
  code_names_["Methionine"]="Methionine";
  code_names_["F"]="Phenylalanine";
  code_names_["PHE"]="Phenylalanine";
  code_names_["Phenylalanine"]="Phenylalanine";
  code_names_["P"]="Proline";
  code_names_["PRO"]="Proline";
  code_names_["Proline"]="Proline";
  code_names_["S"]="Serine";
  code_names_["SER"]="Serine";
  code_names_["Serine"]="Serine";
  code_names_["T"]="Threonine";
  code_names_["THR"]="Threonine";
  code_names_["Threonine"]="Threonine";
  code_names_["W"]="Tryptophan";
  code_names_["TRP"]="Tryptophan";
  code_names_["Tryptophan"]="Tryptophan";
  code_names_["Y"]="Tyrosine";
  code_names_["TYR"]="Tyrosine";
  code_names_["Tyrosine"]="Tyrosine";
  code_names_["V"]="Valine";
  code_names_["VAL"]="Valine";
  code_names_["Valine"]="Valine";
#endif

  //store database infos for copy-ctor
  db_username_ = db_username;
  db_password_ = db_password;
  db_host_     = db_host;

  //fills member variables
  initialize_(type);
}


AminoAcid::~AminoAcid()
{
#ifndef ANNOTATE_XML
  delete sql_adapter_;
#endif
}


AminoAcid::AminoAcid(const AminoAcid& amino_acid)
{
  db_username_= amino_acid.db_username_;
  db_password_= amino_acid.db_password_;
  db_host_    = amino_acid.db_host_;

  //fills member variables
  initialize_(amino_acid.getName());
}


const AminoAcid& AminoAcid::operator=(const AminoAcid& amino_acid)
{
  if (this != &amino_acid) //! no self-assignment
    {
      db_username_= amino_acid.db_username_;
      db_password_= amino_acid.db_password_;
      db_host_    = amino_acid.db_host_;

      //fills member variables
      initialize_(amino_acid.getName());
    }
  return *this;
}


string AminoAcid::getFormula(int position) throw(UnknownPosition)
{
  //position: 0: middle, 1: C-terminal, 2: N-terminal, 3: single

  string pos;
  switch (position)
    {
    case 0:
      return middle_formula_;
    case 1:
      return n_term_formula_;
    case 2:
      return c_term_formula_;
    case 3:
      return single_formula_;
    default:
      throw UnknownPosition(__FILE__, __LINE__, __PRETTY_FUNCTION__,position);
    }
}


double AminoAcid::getMonoMass(int position) throw(UnknownPosition)
{
  //position: 0: middle, 1: N-terminal, 2: C-terminal, 3: single

  string pos;
  switch (position)
    {
    case 0:
      return middle_mono_mass_;
    case 1:
      return n_term_mono_mass_;
    case 2:
      return c_term_mono_mass_;
    case 3:
      return single_mono_mass_;
    default:
      throw UnknownPosition(__FILE__, __LINE__, __PRETTY_FUNCTION__,position);
    }
}


double AminoAcid::getAverageMass(int position) throw(UnknownPosition)
{
  //position: 0: middle, 1: N-terminal, 2: C-terminal, 3: single

  string pos;
  switch (position)
    {
    case 0:
      return middle_average_mass_;
    case 1:
      return n_term_average_mass_;
    case 2:
      return c_term_average_mass_;
    case 3:
      return single_average_mass_;
    default:
      throw UnknownPosition(__FILE__, __LINE__, __PRETTY_FUNCTION__,position);
    }
}


string AminoAcid::getName() const
  {
    return name_;
  }


string AminoAcid::getOneLetter() const
  {
    return one_letter_code_;
  }


string AminoAcid::getThreeLetter() const
  {
    return three_letter_code_;
  }



AminoAcid::UnknownPosition::UnknownPosition(const char* file, int line, const char* function, int request)
throw()
    :	Base(file, line, function, "UnknownPosition", "Requested AminoAcid position not known.")
{
  what_ = ("Position ID \"" + String(request) + "\" not known. (0: middle, 1: C-terminal, 2: N-terminal, 3: single)").c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


AminoAcid::UnknownPosition::~UnknownPosition() throw()
{}


AminoAcid::UnknownAminoAcid::UnknownAminoAcid(const char* file, int line, const char* function, string request)
throw()
    :	Base(file, line, function, "UnknownAminoAcid", "Requested AminoAcid type not known.")
{
  what_ = ("AminoAcid \"" + request + "\" not known.").c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


AminoAcid::UnknownAminoAcid::~UnknownAminoAcid() throw()
{}
