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

#include "Enzyme.h"

#include "ProtDigMembers.h"
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/FORMAT/Param.h>


using namespace std;
using namespace __gnu_cxx;
using namespace OpenMS;


void Enzyme::getID_(string type) throw()
{
  //execute MySQL-query
  sql_adapter_->executeQuery("SELECT enzyme_ID FROM " + (string)ENZ_TABLE + " WHERE enzyme_name  = \"" + type + "\"");

  string ID_string;
  if(sql_adapter_->ifGetUnaryResult(ID_string))
    {
      ID = atoi(ID_string.c_str());
    }
  else
    {
      throw UnknownEnzyme(__FILE__, __LINE__, __PRETTY_FUNCTION__,type);
    }
}


//! constructor
Enzyme::Enzyme()
{
  cerr << "If you use the non-detailed constructor of class Enzyme, you should know what you are doing!" << endl;
}


//! constructor with arguments: \c type specifies what enzyme is used in actual sample
Enzyme::Enzyme(const std::string& ty, string db_username, string db_password, string db_host)
{
  type = ty;

  //store database infos for copy-ctor
  db_username_ = db_username;
  db_password_ = db_password;
  db_host_     = db_host;

#ifndef ANNOTATE_XML
  //connect to database
  sql_adapter_ = new MySQLAdapter();
  sql_adapter_->connect(db_username.c_str(), db_password.c_str(), db_host.c_str());
  sql_adapter_->selectDB(DATABASE);

  //get unique ID from for this instance of AminoAcid from database
  getID_(type);

  //! reads cleavage sites out of database
  sql_adapter_->executeQuery("SELECT cleavage_sites FROM " + (string)ENZ_TABLE + " WHERE enzyme_ID = \"" + String(ID) + "\"");
  string temp_r = sql_adapter_->getUnaryResult();
  for(string::iterator it = temp_r.begin(); it != temp_r.end(); it++)
    {
      cleavage_sites.push_back(string(1,*it)); //! create a string of one copy of char *it
    }

  //! reads reaction mode out of database
  sql_adapter_->executeQuery("SELECT terminality FROM " + (string)ENZ_TABLE + " WHERE enzyme_ID = \"" + String(ID) + "\"");
  cleavage_mode = sql_adapter_->getUnaryResult();
#endif

#ifdef ANNOTATE_XML

  Param param;
  param.load( XML_FILE );
  string temp_r = param.getValue("Preferences:SpecAnnotate:Enzyme:" + type + ":cleav_sites");
  for(string::iterator it = temp_r.begin(); it != temp_r.end(); it++)
    {
      cleavage_sites.push_back(string(1,*it)); //! create a string of one copy of char *it
    }
  cleavage_mode = (string)param.getValue("Preferences:SpecAnnotate:Enzyme:" + type + ":terminality");
#endif

}


//! destructor
Enzyme::~Enzyme()
{}


Enzyme::Enzyme(const Enzyme& enzyme)
{
  type = enzyme.type;

  db_username_= enzyme.db_username_;
  db_password_= enzyme.db_password_;
  db_host_    = enzyme.db_host_;

#ifndef ANNOTATE_XML
  //connect to database
  sql_adapter_ = new MySQLAdapter();
  sql_adapter_->connect(db_username_.c_str(), db_password_.c_str(), db_host_.c_str());
  sql_adapter_->selectDB(DATABASE);

  //get unique ID from for this instance of AminoAcid from database
  getID_(type);

  //! reads cleavage sites out of database
  sql_adapter_->executeQuery("SELECT cleavage_sites FROM " + (string)ENZ_TABLE + " WHERE enzyme_ID = \"" + String(ID) + "\"");
  string temp_r = sql_adapter_->getUnaryResult();
  for(string::iterator it = temp_r.begin(); it != temp_r.end(); it++)
    {
      cleavage_sites.push_back(string(1,*it)); //! create a string of one copy of char *it
    }

  //! reads reaction mode out of database
  sql_adapter_->executeQuery("SELECT terminality FROM " + (string)ENZ_TABLE + " WHERE enzyme_ID = \"" + String(ID) + "\"");
  cleavage_mode = sql_adapter_->getUnaryResult();
#endif

#ifdef ANNOTATE_XML

  Param param;
  param.load( XML_FILE );
  string temp_r = param.getValue("Preferences:SpecAnnotate:Enzyme:" + type + ":cleav_sites");
  for(string::iterator it = temp_r.begin(); it != temp_r.end(); it++)
    {
      cleavage_sites.push_back(string(1,*it)); //! create a string of one copy of char *it
    }
  cleavage_mode = (string) param.getValue("Preferences:SpecAnnotate:Enzyme:" + type + ":terminality");
#endif

}



/* gets necessary pointers to members of \c ProteinDigest, saves fragment start and stop indices in \c ProteinDigest::fragments
  saves to each position all indices of \c ProteinDigest::fragments that signify fragments containing actual position in 
  \c ProteinDigest::sequnce_fragments
 */
void Enzyme::digest(ProtDigMembers members) throw(UnknownCleavageMode)
{
  //! clear sites of application of possible previous instance of \c Enzyme
  members.frags->clear();
  members.cleav_positions->clear();

  //! set an offset concerning cleavage mode of this enzyme
  int offset = 0;
  if (cleavage_mode == "C")
    {
      offset = 1;
    }
  else if (cleavage_mode == "N")
    {
      offset = 0;
    }
  else
    {
      throw UnknownCleavageMode(__FILE__, __LINE__, __PRETTY_FUNCTION__,cleavage_mode);
    }

  /* fill (sortable) list \c ProteinDigest::cleavage_positions:
    iterate over all sites at which this enzyme cleaves
   */
  for(vector<string>::iterator it = cleavage_sites.begin(); it != cleavage_sites.end(); it++)
    {
      AminoAcid temp(*it, db_username_, db_password_, db_host_);
      string threeletter = temp.getThreeLetter(); //! get three letter code of reaction site
      //! does the residue in question occur in sequence of this class
      if (members.aa_occurring->find(threeletter) == members.aa_occurring->end())
        {
          continue; //! no cleavage, since cleavage_site does not occur
        }
      else
        {
          vector<int>::iterator intit;

          //! add all positions of this particular cleavage site to \c ProteinDigest::cleavage_positions
          for(intit = (*(members.aa_positions))[threeletter]->begin(); intit != (*(members.aa_positions))[threeletter]->end(); intit++)
            {
              members.cleav_positions->push_back(*intit);
            }
        }
    }

  //! sort cleavage positions
  members.cleav_positions->sort();

  //! the whole protein is always part of \c fragments
  members.frags->push_back(pair<int,int>(0, (members.seq_oneletter->size()-1)));

  //! add all fragments, that start with beginning of protein:
  for(list<int>
      ::iterator intit = members.cleav_positions->begin();
      intit != members.cleav_positions->end();
      intit++)
    {
      members.frags->push_back(pair<int,int>(0, ((*intit) + (offset-1))));
    }

  //! add all other fragments (starting somewhere in the middle of the protein)
  for(list<int>
      ::iterator a_it = members.cleav_positions->begin();
      a_it != members.cleav_positions->end();
      a_it++)
    {
      list<int>::iterator b_it = a_it;
      for(b_it++; b_it != members.cleav_positions->end(); b_it++)
        {
          members.frags->push_back(pair<int,int>(((*a_it) + offset), ((*b_it) + (offset-1))));
        }
      //! add fragment starting from a_it, going to end of protein
      members.frags->push_back(pair<int,int>(((*a_it) + offset), (members.seq_oneletter->size()-1)));
    }

  //! fill \c ProteinDigest::sequence_fragments
  int x = 0;
  for(vector<pair<int,int> >
      ::iterator it = members.frags->begin();
      it != members.frags->end();
      it++)
    {
      for (int i = it->first; i <= it->second; i++)
        {
          (*members.seq_fragments)[i].push_back(x);
        }
      x++;
    }
}



Enzyme::UnknownEnzyme::UnknownEnzyme(const char* file, int line, const char* function, string request)
throw()
    :	Base(file, line, function, "UnknownEnzyme", "Requested Enzyme type not known.")
{
  what_ = ("Enzyme \"" + request + "\" not known.").c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


Enzyme::UnknownEnzyme::~UnknownEnzyme() throw()
{}


Enzyme::UnknownCleavageMode::UnknownCleavageMode(const char* file, int line, const char* function, string request)
throw()
    :	Base(file, line, function, "UnknownCleavageMode", "Requested cleavage mode not known.")
{
  what_ = ("Cleavage mode \"" + request + "\" not known.").c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


Enzyme::UnknownCleavageMode::~UnknownCleavageMode() throw()
{}
