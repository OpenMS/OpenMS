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
//

#ifndef OPENMS_FORMAT_DATASETINFO_H
#define OPENMS_FORMAT_DATASETINFO_H

#include <OpenMS/FORMAT/PersistentObject.h>

#include <string>
#include <vector>
#include <map>

namespace OpenMS
{
  /**
	  a Dataset is a set of objects in the db<br>
	  a Dataset can include other Datasets<br>
	  circular inclusions are ignored but may lead to errors<br>
	  DataSetInfo include information about Datasets<br>
	  
	  @todo remove
  */
  class DataSetInfo : public PersistentObject
  {
  public:
    /** 
    	@name constructors, destructor, assignment operator <br> 
    */
    //@{
    DataSetInfo();
    DataSetInfo(const DataSetInfo& source);
    DataSetInfo& operator=(const DataSetInfo& source);
    ~DataSetInfo();
    //@}

    int size(std::string type = "") const;
    int dataset_id()const{ return dataset_;}
    const std::vector<int>& contents(std::string type ) const;
    std::string info() const {return info_;}
    std::string name() const {return name_;}

		// Docu in base class
		virtual void persistentWrite(PersistenceManager& pm, const char* name=0) const throw (Exception::Base);
		
		// Docu in base class
		virtual void persistentRead(PersistenceManager& pm) throw (Exception::Base);

	protected:
		// Docu in base class
    virtual void clearChildIds_()
    {
    	//TODO Persistence	
    };		

  private:
    std::string name_;
    std::string info_;
    std::map<std::string,std::vector<int> > contents_;
    int dataset_;
  };
}
#endif //OPENMS_FORMAT_DATASETINFO_H
