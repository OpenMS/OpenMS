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
// $Maintainer:  $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DataSetInfo.h>
#include <OpenMS/FORMAT/PersistenceManager.h>

using namespace std;

namespace OpenMS
{
  DataSetInfo::DataSetInfo()
    :PersistentObject(), name_(),info_(),contents_(),dataset_()
  {
  }

  DataSetInfo::DataSetInfo(const DataSetInfo& source)
    : PersistentObject(source), name_(source.name_),info_(source.info_), contents_(source.contents_),dataset_(source.dataset_)
  {
  }

  DataSetInfo& DataSetInfo::operator=(const DataSetInfo& source)
  {
    PersistentObject::operator=(source);
    name_ = source.name_;
    info_ = source.info_;
    contents_ = source.contents_;
    dataset_ = source.dataset_;
    return *this;
  }
  
  DataSetInfo::~DataSetInfo()
  {
  }

  /**
  \param type
  \return count objects of <i>type</i>, all if type empty
  */
  int DataSetInfo::size(string table) const
  {
    int res = 0;
    if (table == "")
    {
      for (map<string,vector<int> >::const_iterator it = contents_.begin(); it != contents_.end(); ++it)
      {
        res += it->second.size();
      }
    }
    else if ( contents_.find(table) != contents_.end() )
    {
      res = contents_.find(table)->second.size();
    }
    else 
    {
      res = 0;
    }
    return res;
  }

  const vector<int>& DataSetInfo::contents(string type) const
  {
    if ( contents_.find(type) == contents_.end() )
    {
      ostringstream ss;
      ss << "id = " << persistence_id_ << ", type = " << type;
      throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"type not in DataSet",ss.str().c_str());
    }
    return contents_.find(type)->second;
  }

	void DataSetInfo::persistentWrite(PersistenceManager& pm, const char* name) const throw (Exception::Base)
	{
		pm.writeObjectHeader(this,name);
		//TODO Persistence
		pm.writeObjectTrailer(name);
	}
	
	void DataSetInfo::persistentRead(PersistenceManager& pm) throw (Exception::Base)
	{
		//TODO Persistence
		int dummy;
		pm.readPrimitive(dummy,"dummy_");
	}

}
