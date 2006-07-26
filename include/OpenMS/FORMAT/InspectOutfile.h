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
// $Id: InspectOutfile.h,v 1.0 2006/07/12 15:58:59 martinlangwisch Exp $
// $Author: martinlangwisch $
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_INSPECTOUTFILE_H
#define OPENMS_FORMAT_INSPECTOUTFILE_H

#include <OpenMS/FORMAT/Outfile.h>
#include <OpenMS/FORMAT/InspectFile.h>

namespace OpenMS
{
  /**
    @brief Representation of an Inspect outfile
    
    This class serves to read in an Inspect outfile. The information can be 
    retrieved via the >> operator. 
  
  	@ingroup FileIO
  */
  class InspectOutfile:
		public Outfile, 
		InspectFile
  {
    public:
      /// Constructor
      InspectOutfile(const std::string& result_filename, const std::string& database_filename, const std::string& database_path, const double& p_value_threshold = 1.0, std::string index_filename = "") throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument);
      
   };
	
} //namespace OpenMS

#endif // OPENMS_FORMAT_INSPECTOUTFILE_H
