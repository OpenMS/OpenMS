// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_METADATA_ACQUISITION_H
#define OPENMS_METADATA_ACQUISITION_H

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Information about one raw data spectrum that was combined with several
		other raw data spectra.
		
		Although this class is basically a string value, it is needed to store important meta info for each raw data scan. 
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI Acquisition
  	: public MetaInfoInterface
  {
    public:
    	/// Constructor
      Acquisition();
      /// Copy constructor
      Acquisition(const Acquisition& source);
      /// Destructor
      ~Acquisition();
      
      /// Assignment operator
      Acquisition& operator= (const Acquisition& source);
			
			/// Equality operator
      bool operator== (const Acquisition& rhs) const;
      /// Equality operator
      bool operator!= (const Acquisition& rhs) const;
			
			/// return the identifier/index/number of the acquisition
      const String& getIdentifier() const;
      /// sets the index/number of the scan
      void setIdentifier(const String& identifier);

    protected:
      String identifier_;
      
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_ACQUISITION_H
