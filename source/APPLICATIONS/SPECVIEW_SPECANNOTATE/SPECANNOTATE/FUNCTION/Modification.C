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



#include "Modification.h"
#include "ProtDigMembers.h"

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/Param.h>


using namespace std;



void Modification::getID_(string type) throw()
{
  //execute MySQL-query
  sql_adapter_->executeQuery("SELECT modification_ID FROM " + (string)MOD_TABLE + " WHERE modification_name  = \"" + type + "\"");

  string ID_string;
  if(sql_adapter_->ifGetUnaryResult(ID_string))
    {
      ID = atoi(ID_string.c_str());
    }
  else
    {
      throw UnknownModification(__FILE__, __LINE__, __PRETTY_FUNCTION__,type);
    }
}


void Modification::initialize_()
{
#ifdef ANNOTATE_XML
  Param param;
  param.load( XML_FILE );
  plus_formula = (string)param.getValue("Preferences:SpecAnnotate:Modification:" + String(ID) + ":plus_formula");
  plus_mono_mass = (string)param.getValue("Preferences:SpecAnnotate:Modification:" + String(ID) + ":plus_mono_mass");
  minus_mono_mass = (string)param.getValue("Preferences:SpecAnnotate:Modification:" + String(ID) + ":minus_mono_mass");
  plus_average_mass = (string)param.getValue("Preferences:SpecAnnotate:Modification:" + String(ID) + ":plus_average_mass");
  minus_average_mass = (string)param.getValue("Preferences:SpecAnnotate:Modification:" + String(ID) + ":minus_average_mass");
#endif

#ifndef ANNOTATE_XML
  //! reads reaction sites out of inifile
  sql_adapter_->executeQuery("SELECT modification_sites FROM "
                             + (string)MOD_TABLE + " WHERE modification_ID = \"" + String(ID) + "\"");
  string temp_r = sql_adapter_->getUnaryResult();
  for(string::iterator it = temp_r.begin(); it != temp_r.end(); it++)
    {
      modification_sites.push_back(string(1,*it)); //! create a string of one copy of char *it
    }

  //! fill member variables
  sql_adapter_->executeQuery("SELECT plus_formula FROM "
                             + (string)MOD_TABLE + " WHERE modification_ID = \"" + String(ID) + "\"");
  plus_formula = sql_adapter_->getUnaryResult();

  sql_adapter_->executeQuery("SELECT minus_formula FROM "
                             + (string)MOD_TABLE + " WHERE modification_ID = \"" + String(ID) + "\"");
  minus_formula = sql_adapter_->getUnaryResult();

  sql_adapter_->executeQuery("SELECT plus_mono_mass FROM "
                             + (string)MOD_TABLE + " WHERE modification_ID = \"" + String(ID) + "\"");
  plus_mono_mass = sql_adapter_->getUnaryResult();

  sql_adapter_->executeQuery("SELECT minus_mono_mass FROM "
                             + (string)MOD_TABLE + " WHERE modification_ID = \"" + String(ID) + "\"");
  minus_mono_mass = sql_adapter_->getUnaryResult();

  sql_adapter_->executeQuery("SELECT plus_average_mass FROM "
                             + (string)MOD_TABLE + " WHERE modification_ID = \"" + String(ID) + "\"");
  plus_average_mass = sql_adapter_->getUnaryResult();

  sql_adapter_->executeQuery("SELECT minus_average_mass FROM "
                             + (string)MOD_TABLE + " WHERE modification_ID = \"" + String(ID) + "\"");
  minus_average_mass = sql_adapter_->getUnaryResult();
#endif
}


//! constructor
Modification::Modification()
{
  cerr << "If you use the non-detailed constructor of class Modification, you should know what you are doing!" << endl;
}


//! constructor with argument: \c id specifies the database ID of the modification (if already known)
Modification::Modification(int id, string db_username, string db_password, string db_host)
{
  ID = id;

  //store database infos for copy-ctor
  db_username_ = db_username;
  db_password_ = db_password;
  db_host_     = db_host;

#ifdef ANNOTATE_XML

  Param param;
  param.load( XML_FILE );
  type = (string)param.getValue("Preferences:SpecAnnotate:Modification:" + String(id) + ":name");
#endif

#ifndef ANNOTATE_XML
  //connect to database
  sql_adapter_ = new MySQLAdapter();
  sql_adapter_->connect(db_username.c_str(), db_password.c_str(), db_host.c_str());
  sql_adapter_->selectDB(DATABASE);

  //get type
  sql_adapter_->executeQuery("SELECT modification_name FROM " + (string)MOD_TABLE + " WHERE modification_ID  = \""
                             + String(id) + "\"");
  type = sql_adapter_->getUnaryResult();
#endif

  //read member variables
  initialize_();
}


//! constructor with argument: \c type specifies the type of the modification
Modification::Modification(const std::string& ty, string db_username, string db_password, string db_host)
{

#ifdef ANNOTATE_XML
  cerr << "With no database present, class Modification only can be initialized via ID, not type/name" << endl;
  throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong initialization of class Modification", "No ID specified!");
#endif

  type = ty;

  //store database infos for copy-ctor
  db_username_ = db_username;
  db_password_ = db_password;
  db_host_     = db_host;

  //connect to database
  sql_adapter_ = new MySQLAdapter();
  sql_adapter_->connect(db_username.c_str(), db_password.c_str(), db_host.c_str());
  sql_adapter_->selectDB(DATABASE);

  //get unique ID from for this instance of Modification from database
  getID_(type);

  //read member variables
  initialize_();
}


//! copy constructor: explicitly initialize base class, an appropriate assignment operator ist generated automatically
Modification::Modification(const Modification& mod)
{
  // for new instance of \c Modification no database connection is needed, alle information already read from database

  type = mod.type;

  ID = mod.ID;

  for(vector<string>::const_iterator it = (mod.modification_sites).begin(); it != (mod.modification_sites).end(); it++)
    {
      modification_sites.push_back(*it);
    }

  plus_formula = mod.plus_formula;
  minus_formula = mod.minus_formula;

  plus_mono_mass = mod.plus_mono_mass;
  minus_mono_mass = mod.minus_mono_mass;

  plus_average_mass = mod.plus_average_mass;
  minus_average_mass = mod.minus_average_mass;
}


//! modifies amino amino acids specified in \c modification_sites . gets pointers to necessarry members of \c ProteinDigest
void Modification::modifyOverall(ProtDigMembers members)
{
  //! iterate over all sites to which this modification applies
  for(vector<string>::iterator it = modification_sites.begin(); it != modification_sites.end(); it++)
    {
      AminoAcid temp(*it, db_username_, db_password_, db_host_);
      string threeletter = temp.getThreeLetter(); //! get three letter code of reaction site
      //! does the residue in question occur in sequence of this class
      if (members.aa_occurring->find(threeletter) != members.aa_occurring->end())
        {
          vector<int>::iterator intit;
          for(intit = (*(members.aa_positions))[threeletter]->begin(); intit != (*(members.aa_positions))[threeletter]->end(); intit++)
            {
              if ((*(members.seq_overall_modifications))[(*intit)] == 0)
                {
                  (*(members.seq_overall_modifications))[(*intit)] = ID;
                }
              else
                {
                  throw AmbiguousOverallModification(__FILE__, __LINE__, __PRETTY_FUNCTION__,type);
                }
            }
        }
    }
}


bool Modification::canModify(string residue)
{
  if (find(modification_sites.begin(), modification_sites.end(), residue) != modification_sites.end())
    {
      return true;
    }
  else
    {
      return false;
    }
}



void Modification::changeID(int new_ID)
{
  ID = new_ID;

#ifndef ANNOTATE_XML
  //get type
  sql_adapter_->executeQuery("SELECT modification_name FROM " + (string)MOD_TABLE + " WHERE modification_ID  = \""
                             + String(ID) + "\"");
  type = sql_adapter_->getUnaryResult();
#endif

#ifdef ANNOTATE_XML

  Param param;
  param.load( XML_FILE );
  type = (string)param.getValue("Preferences:SpecAnnotate:Modification:" + String(ID) + ":name");
#endif

  //read member variables
  initialize_();
}


//! just relaying to \c MolecularFormula::getMonoMass(), \c sign: 0: + formula, 1 - formula
float Modification::getMonoMass(int sign) throw()
{
  if (sign == 0)
    {
      return atof(plus_mono_mass.c_str());
    }
  else if (sign == 1)
    {
      return atof(minus_mono_mass.c_str());
    }
  else
    {
      throw UnknownFormula(__FILE__, __LINE__, __PRETTY_FUNCTION__,sign);
    }
}


float Modification::getAverageMass(int sign) throw()
{
  if (sign == 0)
    {
      return atof(plus_average_mass.c_str());
    }
  else if (sign == 1)
    {
      return atof(minus_average_mass.c_str());
    }
  else
    {
      throw UnknownFormula(__FILE__, __LINE__, __PRETTY_FUNCTION__,sign);
    }
}


string Modification::getMolecularFormula(int sign) throw()
{
  if (sign == 0)
    {
      return plus_formula;
    }
  else if (sign == 1)
    {
      return minus_formula;
    }
  else
    {
      throw UnknownFormula(__FILE__, __LINE__, __PRETTY_FUNCTION__,sign);
    }
}


Modification::UnknownModification::UnknownModification(const char* file, int line, const char* function, string request)
throw()
    :	Base(file, line, function, function, "UnknownModification", "Requested Modification type not known.")
{
  what_ = ("Modification \"" + request + "\" not known.").c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


Modification::UnknownModification::~UnknownModification() throw()
{}


Modification::UnknownFormula::UnknownFormula(const char* file, int line, const char* function, int request)
throw()
    :	Base(file, line, function, function, "UnknownFormula", "Requested Formula not known.")
{
  what_ = ("Requested Formula, accession code:  \"" + String(request) + "\", not known. (0: + formula, 1: - formula)").c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


Modification::UnknownFormula::~UnknownFormula() throw()
{}


Modification::AmbiguousOverallModification::AmbiguousOverallModification(const char* file, int line, const char* function, std::string type)
throw()
    :	Base(file, line, function, function, "AmbiguousOverallModification", "Two modifications claim same residue.")
{
  what_ = ("Overall Modification \"" + type + "\" claims residue, that is already modified").c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


Modification::AmbiguousOverallModification::~AmbiguousOverallModification() throw()
{}
