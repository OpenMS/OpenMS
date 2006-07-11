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


#include "ModificationStringParser.h"




//namespaces
using namespace std;
using namespace OpenMS;


//! default constructor
ModificationStringParser::ModificationStringParser(std::string db_username, std::string db_password, std::string db_host)
{
  db_username_ = db_username;
  db_password_ = db_password;
  db_host_ = db_host;
}


bool ModificationStringParser::readPosition(istringstream& temp_string_stream, int& position)
{
  string tmp_string;
  temp_string_stream >> tmp_string;
  position = atoi(tmp_string.c_str());

  //test if position valid
  if ((tmp_string != "0") && (position == 0))
    {
      return false;
    }
  return true;
}


//! reads number of occurrence out of \c temp_string
bool ModificationStringParser::readModification(istringstream& temp_string_stream, int& mod_ID)
{
  string tmp_string;
  temp_string_stream >> tmp_string;
  mod_ID = atoi(tmp_string.c_str());

  //test if mod_ID valid
  if (mod_ID == 0)
    {
      return false;
    }
  return true;
}


//! does the work
vector<pair<int, vector<Modification*> > > ModificationStringParser::parse(const string& mod_string) throw(InvalidModificationString)
{
  vector<pair<int, vector<Modification*> > > result;
  int position, mod_ID;
  Modification* mod;
  vector<Modification*> mods;
  istringstream temp_string_stream(mod_string);
  string tmp;

  for (;!temp_string_stream.eof();)
    {
      //! try to read a position out of temp_string_stream
      if ( !readPosition(temp_string_stream, position) )
        {
          throw InvalidModificationString(__FILE__, __LINE__, __PRETTY_FUNCTION__,mod_string);
        }

      //! try to read set of modifications out of temp_string_stream
      temp_string_stream >> tmp;
      if (tmp != "(")
        {
          throw InvalidModificationString(__FILE__, __LINE__, __PRETTY_FUNCTION__,mod_string);
        }
      for (;;)
        {
          if (temp_string_stream.eof())
            {
              throw InvalidModificationString(__FILE__, __LINE__, __PRETTY_FUNCTION__,mod_string);
            }
          if ( !readModification(temp_string_stream, mod_ID) )
            {
              throw InvalidModificationString(__FILE__, __LINE__, __PRETTY_FUNCTION__,mod_string);
            }

          mod = new Modification(mod_ID, db_username_, db_password_, db_host_);
          mods.push_back(mod);

          temp_string_stream >> tmp;
          if (tmp == ")")
            {
              break;
            }
          else if (tmp != ",")
            {
              throw InvalidModificationString(__FILE__, __LINE__, __PRETTY_FUNCTION__,mod_string);
            }
        }

      result.push_back(pair<int, vector<Modification*> >(position, mods));
      mods.clear();

      temp_string_stream >> tmp;
      if (tmp == "*")
        {
          break;
        }
      else if (tmp == ";") //just not to leave a case out
        {
          continue;
        }
    }
  return result;
}



ModificationStringParser::InvalidModificationString::InvalidModificationString(const char* file, int line, const char* function, string mod_string)
throw()
    :	Base(file, line, function, function, "InvalidModificationString", "Given string for partial modifcations invalid.")
{
  what_ = (mod_string + " is not a valid partial modification string!").c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


ModificationStringParser::InvalidModificationString::~InvalidModificationString() throw()
{}
