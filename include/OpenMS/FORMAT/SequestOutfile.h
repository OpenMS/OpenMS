/// -*- Mode: C++; tab-width: 2; -*-
/// vi: set ts=2:
///
/// --------------------------------------------------------------------------
///                   OpenMS Mass Spectrometry Framework
/// --------------------------------------------------------------------------
///  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
///
///  This library is free software; you can redistribute it and/or
///  modify it under the terms of the GNU Lesser General Public
///  License as published by the Free Software Foundation; either
///  version 2.1 of the License, or (at your option) any later version.
///
///  This library is distributed in the hope that it will be useful,
///  but WITHOUT ANY WARRANTY; without even the implied warranty of
///  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
///  Lesser General Public License for more details.
///
///  You should have received a copy of the GNU Lesser General Public
///  License along with this library; if not, write to the Free Software
///  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
///
/// --------------------------------------------------------------------------
/// $Id: SequestOutfile.h,v 1.0 2006/07/12 15:58:59 martinlangwisch Exp $
/// $Author: martinlangwisch $
/// $Maintainer: Martin Langwisch $
/// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_SEQUESTOUTFILE_H
#define OPENMS_FORMAT_SEQUESTOUTFILE_H


/*#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/DATASTRUCTURES/Date.h>*/
#include <OpenMS/FORMAT/Outfile.h>

#include <sstream>
#include <string>
#include <map>
#include <set>
#include <cstring>

namespace OpenMS
{
  /**
    @brief Representation of an Sequest outfile
    
    This class serves to read in an Sequest outfile. The information can be 
    retrieved via the >> operator. 
  
  	@ingroup FileIO
  */
  class SequestOutfile:
		public Outfile
  {
    public:
			/// empty constructor
			SequestOutfile();
			
			/// copy constructor
			SequestOutfile(const SequestOutfile& source);
			
      /// Constructor
      SequestOutfile(const std::string& result_filename) throw (Exception::FileNotFound, Exception::ParseError);
			
		protected:
			bool split_(const String& s, const String& splitter, std::vector<String>& substrings);
			
			void getColumns_(std::ifstream& result_file, std::vector< String >& substrings, String& result_filename);
			
			std::set< String > getSequences_(const String& database_path, const String& database_filename, std::set< String > acdt_set, std::vector< String >& sequences);
   };
	
} //namespace OpenMS

#endif // OPENMS_FORMAT_SEQUESTOUTFILE_H
