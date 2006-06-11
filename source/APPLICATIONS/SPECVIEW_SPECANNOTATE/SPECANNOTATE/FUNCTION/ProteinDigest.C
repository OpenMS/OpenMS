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
// $Id: ProteinDigest.C,v 1.4 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------



#include "ProteinDigest.h"
#include "ProtDigMembers.h"

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/Param.h>

using namespace std;
using namespace __gnu_cxx;



//! default constructor
ProteinDigest::ProteinDigest()
{}


ProteinDigest::ProteinDigest(string identifier, int ID, string db_username, string db_password, string db_host)
{
  //store database infos for copy-ctor
  db_username_ = db_username;
  db_password_ = db_password;
  db_host_     = db_host;

  protein_identifier_ = identifier;
  id = ID;

  //! read and store all necessary information
  initialize();
}


//! copy constructor
ProteinDigest::ProteinDigest(const ProteinDigest& protein_digest)
{
  //store database infos
  db_username_ = protein_digest.db_username_;
  db_password_ = protein_digest.db_password_;
  db_host_     = protein_digest.db_host_;

  protein_identifier_ = protein_digest.protein_identifier_;

  initialize();
}


//! assignment operator
const ProteinDigest& ProteinDigest::operator=(const ProteinDigest& protein_digest)
{
  if (this != &protein_digest)
    {
      //store database infos
      db_username_ = protein_digest.db_username_;
      db_password_ = protein_digest.db_password_;
      db_host_     = protein_digest.db_host_;

      protein_identifier_ = protein_digest.protein_identifier_;

      initialize();
    }
  return *this;
}


/* initializes different representations of sequence of protein to be digested in synchronous way
  fills \c protein , \c sequence_oneletter , \c sequence_threeletter and \c sequence_aminoacids with the right values.
  sets \c sequence_fragments to be of the right size.
  this method futhermore fills \c aminoacids_occurring and \c aminacids_positions with the right values
 */
void ProteinDigest::initialize() throw()
{

#ifndef ANNOTATE_XML
  //connect to database
  sql_adapter_ = new MySQLAdapter();
  sql_adapter_->connect(db_username_.c_str(), db_password_.c_str(), db_host_.c_str());
  sql_adapter_->selectDB(DATABASE);

  //get ID of Protein in Database
  sql_adapter_->executeQuery("SELECT protein_ID FROM " + (string)PROTEIN_TABLE + " WHERE identifier = \"" + protein_identifier_ + "\"");
  protein_ID_=atoi((sql_adapter_->getUnaryResult()).c_str());

  //try to get sequence of protein out of database
  sql_adapter_->executeQuery("SELECT sequence_oneletter FROM " + (string)PROTEIN_TABLE + " WHERE protein_ID = \""
                             + String(protein_ID_) + "\"");
  sequence_oneletter=sql_adapter_->getUnaryResult();

  //if that`s not possible, get sequence out of file
  if (sequence_oneletter == "")
    {
      //get FASTA-Filename of Protein out of Database
      sql_adapter_->executeQuery("SELECT fasta_filename FROM " + (string)PROTEIN_TABLE + " WHERE protein_ID = \""
                                 + String(protein_ID_) + "\"");
      protein_filename=sql_adapter_->getUnaryResult();

      if (protein_filename == "void")
        {
          throw NoProteinFilename(__FILE__, __LINE__, __PRETTY_FUNCTION__,protein_identifier_);
        }

      //! reading out of FASTA File
      ifstream infile(protein_filename.c_str());
      if (!infile)
        {
          throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong filename", ("Could not open file " + protein_filename).c_str() );
        }

      // throw first line away: not interesting
      string line;
      getline(infile,line);

      //iterate rest of lines: sequence_oneletter
      for (getline(infile,line);!infile.eof();getline(infile,line))
        {
          istringstream ist(line);
          string tmp;
          ist >> tmp;

          sequence_oneletter += tmp;
        }
    }
#endif

#ifdef ANNOTATE_XML
  Param param;
  param.load( XML_FILE );
  sequence_oneletter = (string)param.getValue("Preferences:SpecAnnotate:Protein:" + protein_identifier_ + ":sequence_oneletter");
#endif

  //! generating \c sequence_aminoacids and \c sequence_threeletter at once
  for (string::iterator it = sequence_oneletter.begin(); it != sequence_oneletter.end(); it++)
    {
      AminoAcid* tmp_aa = new AminoAcid((String(*it)),db_username_, db_password_, db_host_);
      string threeletter = tmp_aa->getThreeLetter();

      //! if amino acid didn't occur yet, insert temporary instance of \c AminoAcid
      if (aminoacids_occurring.find(threeletter) == aminoacids_occurring.end())
        {
          aminoacids_occurring[threeletter] = tmp_aa;
          aminoacids_positions[threeletter] = new vector<int>;
        }
      else //clean up
        {
          delete tmp_aa;
        }

      sequence_aminoacids.push_back(aminoacids_occurring[threeletter]);        //! build up \c sequence_aminoacids
      sequence_threeletter.push_back(threeletter); //! build up \c sequence_threeletter

      //! build up \c aminoacid_positions: -1 since positions start with 0
      aminoacids_positions[threeletter]->push_back(sequence_threeletter.size()-1);
    }

  //! set rest of sequence representations respectivley info-containers to be of the right size
  sequence_fragments.resize(sequence_threeletter.size());
  sequence_overall_modifications.resize(sequence_threeletter.size(),0);

  //! be sure instances of following members are created:
  fragments.clear();
  cleavage_positions.clear();


#ifndef ANNOTATE_XML
  //fill some information into database values are calculated taking given filename into account. already present values will be overwritten
  //insert sequence_oneletter into database
  sql_adapter_->executeQuery("UPDATE " + (string)PROTEIN_TABLE + " SET `sequence_oneletter` = \"" + sequence_oneletter + "\" WHERE "
                             + " `protein_ID` = \"" + String(protein_ID_) + "\"");

  //insert number of aminoacids
  sql_adapter_->executeQuery("UPDATE " + (string)PROTEIN_TABLE + " SET `no_of_aminoacids` = \"" + String(sequence_oneletter.size())
                             + "\" WHERE `protein_ID` = \"" + String(protein_ID_) + "\"");

  //insert masses
  sql_adapter_->executeQuery("UPDATE " + (string)PROTEIN_TABLE + " SET `mono_mass` = \""
                             + String(getFragmentMonoMass(0, (sequence_oneletter.size()-1)))
                             + "\" WHERE `protein_ID` = \"" + String(protein_ID_) + "\"");
  sql_adapter_->executeQuery("UPDATE " + (string)PROTEIN_TABLE + " SET `average_mass` = \""
                             + String(getFragmentAverageMass(0, (sequence_oneletter.size()-1)))
                             + "\" WHERE `protein_ID` = \"" + String(protein_ID_) + "\"");

  //fill table with sequences
  int pos = 0;
  for(vector<string>::iterator it = sequence_threeletter.begin(); it != sequence_threeletter.end(); it++)
    {
      //get aminoacid id for this position
      sql_adapter_->executeQuery("SELECT aminoacid_ID FROM " + (string)AMINO_TABLE + " WHERE three_letter_code = \"" + (*it) + "\"");
      string amino_id = sql_adapter_->getUnaryResult();

      //if entry doesn't exist yet, insert it
      sql_adapter_->executeQuery("SELECT count(*) FROM " + (string)SEQUENCE_TABLE + " WHERE protein_ID = " + String(protein_ID_)
                                 + " AND s_position = " + String(pos));
      if (atoi((sql_adapter_->getUnaryResult()).c_str()) == 0)
        {
          //insert
          sql_adapter_->executeQuery("INSERT INTO " + (string)SEQUENCE_TABLE + " ( `protein_ID` , `s_position` , `aminoacid_ID` ) "
                                     + " VALUES ( '" + String(protein_ID_) + "', '" + String(pos) + "', '" + amino_id +"' )");
        }

      //go to next position
      pos++;
    }
#endif
}



//! via this method instances of \c Enzyme can be applied to contents of this class
void ProteinDigest::digest(Enzyme* enz)
{
  //! create class with necessary members for \c Enzyme
  ProtDigMembers members(this);

  //! digest!
  enz->digest(members);

  //! fill table digest_fragment of database and retrieve database index of each fragment (if already present there)
  for (vector<pair<int,int> >::iterator it = fragments.begin(); it != fragments.end(); it++)
    {
      //if entry doesn't exist yet, insert it
      sql_adapter_->executeQuery("SELECT count(*) FROM " + (string)FRAGMENT_TABLE + " WHERE protein_ID = " + String(protein_ID_)
                                 + " AND enzyme_ID   = " + String(enz->getID())
                                 + " AND d_start_pos = " + String(it->first)
                                 + " AND d_end_pos   = " + String(it->second));
      if (atoi((sql_adapter_->getUnaryResult()).c_str()) == 0)
        {
          //insert
          sql_adapter_->executeQuery("INSERT INTO " + (string)FRAGMENT_TABLE
                                     + " ( `protein_ID` , `enzyme_ID` , `d_start_pos` , `d_end_pos` ) "
                                     + " VALUES ( '" + String(protein_ID_) + "', '" + String(enz->getID())
                                     + "', '" + String(it->first) + "', '" + String(it->second) + "' )");
        }

      //add database index of fragment to \c fragment_database_indices_
      sql_adapter_->executeQuery("SELECT digest_fragment_ID FROM " + (string)FRAGMENT_TABLE + " WHERE protein_ID = "
                                 + String(protein_ID_)
                                 + " AND enzyme_ID   = " + String(enz->getID())
                                 + " AND d_start_pos = " + String(it->first)
                                 + " AND d_end_pos   = " + String(it->second));
      fragment_database_indices_.push_back(atoi((sql_adapter_->getUnaryResult()).c_str()));

    }
}



void ProteinDigest::modifyOverall(Modification* mod)
{
  //! create class with necessary members for \c Enzyme
  ProtDigMembers members(this);

  //! modify overall
  mod->modifyOverall(members);
}


int ProteinDigest::dbStoreOverallModifications()
{
  // add set of overall modifiications to database (table realized_modifications)
  string first_mod_ID = "-1";       //for returning database ID of first modification
  string last_mod_ID = "-1";        //for getting pointers right
  string actual_mod_ID;
  for(unsigned int i=1; i < sequence_overall_modifications.size(); i++)
    {
      // residue modified?
      if (sequence_overall_modifications[i] != 0)
        {
          // add actual position/modification
          sql_adapter_->executeQuery("INSERT INTO " + (string)REALIZED_MOD_TABLE
                                     + " ( `m_position` , `modification_ID` ) "
                                     + " VALUES ( '" + String(i) + "', '" + String(sequence_overall_modifications[i]) + "' )");

          //get database id of actual modification
          sql_adapter_->executeQuery("SELECT last_insert_id() FROM " + (string)REALIZED_MOD_TABLE + " LIMIT 1");
          actual_mod_ID = sql_adapter_->getUnaryResult();


          // if we do not work on the first modification: add pointer to actual modification into last modification
          if (last_mod_ID != "-1")
            {
              sql_adapter_->executeQuery("UPDATE " + (string)REALIZED_MOD_TABLE
                                         + " SET `next_realized_modification_ID` = " + actual_mod_ID
                                         + " WHERE `realized_modification_ID` = " + last_mod_ID + " LIMIT 1");
            }


          //store this position
          last_mod_ID = actual_mod_ID;

          //if this is first modification: store in variable.
          if (first_mod_ID == "-1")
            {
              first_mod_ID = actual_mod_ID;
            }
        }
    }
  return atoi(first_mod_ID.c_str());
}



/* returns list of vectors of ints:
  first elements of vectors are indices in \c fragments of that fragments, that contain certain a site (given in one-letter-code)
  second elements of vectors are positions of sites
  third elements of vectors are ID's of this instance of \c ProteinDigest
  fourth elements of vectors are indices of the fragments in database (may or may not be the same ase in \c fragments)
 */
list<vector<int> > ProteinDigest::getFragmentIndicesContaining(const vector<string>& sites)
{
  list<int> positions;

  //! iterate over sites and fill positions
  for(vector<string>::const_iterator it = sites.begin(); it != sites.end(); it++)
    {
      AminoAcid temp(*it, db_username_, db_password_, db_host_);
      string threeletter = temp.getThreeLetter(); //! get three letter code of link site
      //! does the residue in question occur in sequence of this class
      if (aminoacids_occurring.find(threeletter) != aminoacids_occurring.end())
        {
          vector<int>::iterator intit;
          for(intit = aminoacids_positions[threeletter]->begin(); intit != aminoacids_positions[threeletter]->end(); intit++)
            {
              positions.push_back(*intit);
            }
        }
    }

  //! iterate over all link positions and fill result every entry in \c result: first int: frag index, second int: link pos, third int: id
  list<vector<int> > result;
  for(list<int>::iterator it = positions.begin(); it != positions.end(); it++)
    {
      for(vector<int>::iterator intit = sequence_fragments[*it].begin(); intit != sequence_fragments[*it].end(); intit++)
        {
          vector<int> temp;
          temp.push_back(*intit);                              //! fragment id
          temp.push_back(*it);                                 //! link position
          temp.push_back(id);                                  //! id of this instance of \c ProteinDigest
          temp.push_back(fragment_database_indices_[*intit]);  //! id of this fragment in database (table digest_fragments)
          result.push_back(temp);
        }
    }

  return result;
}



double ProteinDigest::getFragmentOverallModifiedMass(int start_pos, int end_pos, string masstype, hash_map<int,int>& temp_o_mods,
    Modification& mod)
{
  //determine no of occurrences of overall modifications
  for (int i = start_pos; i <= end_pos; i++)
    {
      //same overall modification already seen?
      if (temp_o_mods.find(sequence_overall_modifications[i]) != temp_o_mods.end())
        {
          (temp_o_mods[(sequence_overall_modifications[i])])++;
        }
      //first occurrence of overall mod with id *it
      else if (sequence_overall_modifications[i] != 0)
        {
          temp_o_mods[(sequence_overall_modifications[i])] = 1;
        }
    }

  //get unmodified mass of fragment
  double frag_ov_mod_mass = 0;
  if (masstype == "average")
    {
      frag_ov_mod_mass = getFragmentAverageMass(start_pos, end_pos);
    }
  else if (masstype == "mono")
    {
      frag_ov_mod_mass = getFragmentMonoMass(start_pos, end_pos);
    }
  else
    {
      throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "No proper masstype",
                                    (masstype + " is unknown masstype.").c_str());
    }

  //add masses of overall modifications
  for (hash_map<int,int>::iterator it = temp_o_mods.begin(); it != temp_o_mods.end(); it++)
    {
      mod.changeID(it->first);
      if (masstype == "average")
        {
          frag_ov_mod_mass += ((it->second) * (mod.getAverageMass(0)));
          frag_ov_mod_mass -= ((it->second) * (mod.getAverageMass(1)));
        }
      else if (masstype == "mono")
        {
          frag_ov_mod_mass += ((it->second) * (mod.getMonoMass(0)));
          frag_ov_mod_mass -= ((it->second) * (mod.getMonoMass(1)));
        }
    }
  return frag_ov_mod_mass;
}



//! returns mass of fragment given by start and end position. start and end aa's are considered to be N-/C-terminal (WITHOUT modifications)
double ProteinDigest::getFragmentMonoMass(int start_pos, int end_pos)
{
  double result = 0;
  if (start_pos == end_pos)
    {
      //! a single amino acid (accession code 3)
      result += sequence_aminoacids[start_pos]->getMonoMass(3);
    }
  else if (start_pos <= end_pos)
    {

      //! first position is N-terminal (accession code 1)
      result += sequence_aminoacids[start_pos]->getMonoMass(1);

      //! positions in the middle (accession code 0)
      for (int i = (start_pos+1); i < end_pos; i++)
        {
          result += sequence_aminoacids[i]->getMonoMass(0);
        }

      //! last position is C-terminal (accession code 2)
      result += sequence_aminoacids[end_pos]->getMonoMass(2);

    }
  else
    {
      throw WrongPositionInProtein(__FILE__, __LINE__, __PRETTY_FUNCTION__,"getFragmentMonoMass()");
    }

  return result;
}



//! returns mass of fragment given by start and end position. start and end aa's are considered to be N-/C-terminal
double ProteinDigest::getFragmentAverageMass(int start_pos, int end_pos)
{
  double result = 0;
  if (start_pos == end_pos)
    {
      //! a single amino acid (accession code 3)
      result += sequence_aminoacids[start_pos]->getAverageMass(3);
    }
  else if (start_pos <= end_pos)
    {

      //! first position is N-terminal (accession code 1)
      result += sequence_aminoacids[start_pos]->getAverageMass(1);

      //! positions in the middle (accession code 0)
      for (int i = (start_pos+1); i < end_pos; i++)
        {
          result += sequence_aminoacids[i]->getAverageMass(0);
        }

      //! last position is C-terminal (accession code 2)
      result += sequence_aminoacids[end_pos]->getAverageMass(2);

    }
  else
    {
      throw WrongPositionInProtein(__FILE__, __LINE__, __PRETTY_FUNCTION__,"getFragmentAverageMass()");
    }

  return result;
}



ProteinDigest::NoProteinFilename::NoProteinFilename(const char* file, int line, const char* function, string request)
throw()
    :	Base(file, line, function, "NoProteinFilename", "No fasta-filename present in database.")
{
  what_ = ("For identifier \"" + request + "\" is no fasta-filename present in database.").c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


ProteinDigest::NoProteinFilename::~NoProteinFilename() throw()
{}


ProteinDigest::WrongPositionInProtein::WrongPositionInProtein(const char* file, int line, const char* function, string method)
throw()
    :	Base(file, line, function, "WrongPositionInProtein", "Wrong start- or end-position given.")
{
  what_ = ("Wrong start- or end-position of fragment given in method" + method).c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


ProteinDigest::WrongPositionInProtein::~WrongPositionInProtein() throw()
{}


ProteinDigest::WrongModification::WrongModification(const char* file, int line, const char* function, int mod_ID, int pos)
throw()
    :	Base(file, line, function, "WrongModification", "Given modification not able to modify residue at given position.")
{
  what_ = ("Modification with ID " + String(mod_ID) + " not able to modify residue at position " + String(pos)).c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


ProteinDigest::WrongModification::~WrongModification() throw()
{}
