// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_ACQUISITIONINFO_H
#define OPENMS_METADATA_ACQUISITIONINFO_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Acquisition.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <vector>

namespace OpenMS 
{
	/**
		@brief Description of the combination of raw data to a single spectrum
		
		Specification for combining raw scans ( Acquisition ) into a single spectrum. 
		A list of acquisitions from the original raw file can be specified.
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI AcquisitionInfo
  	: public std::vector<Acquisition>,
  		public MetaInfoInterface
  {
    public:
    	/// Constructor
      AcquisitionInfo();
      /// Copy constructor
      AcquisitionInfo(const AcquisitionInfo& source);
      /// Destructor
      ~AcquisitionInfo();
			
			/// Assignment operator
      AcquisitionInfo& operator= (const AcquisitionInfo& source);

      /// Equality operator
      bool operator== (const AcquisitionInfo& rhs) const;
      /// Equality operator
      bool operator!= (const AcquisitionInfo& rhs) const;
			
			/// returns the method of combination
      const String& getMethodOfCombination() const;
      /// sets the method of combination
      void setMethodOfCombination(const String& method_of_combination);
			
    protected:
      String method_of_combination_;
      
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_ACQUISITIONINFO_H
