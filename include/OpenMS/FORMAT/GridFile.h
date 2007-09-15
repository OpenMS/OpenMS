// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_GRIDFILE_H
#define OPENMS_FORMAT_GRIDFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/Grid.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/GridHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

namespace OpenMS
{

  /**
    @brief Provides Input/Output functionality for instances of class DGrid.
    
    @ingroup FileIO
   */
  class GridFile
		: public Internal::XMLFile
  {
	  public:
	    ///Default constructor
	    GridFile()
	    {
	    }
	
	    ///Destructor
	    virtual ~GridFile()
	    {
	    }
	
	    /// loads the file with name @p filename into @p grid.
	    void load(String filename, Grid& grid) throw (Exception::FileNotFound,Exception::ParseError);
	
	
	    /// stores the grid @p grid in file with name @p filename.
	    void store(String filename, const Grid& grid) const throw (Exception::UnableToCreateFile);
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_DGRIDFILE_H
