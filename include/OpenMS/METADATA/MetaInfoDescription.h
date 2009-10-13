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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_METAINFODESCRIPTION_H
#define OPENMS_METADATA_METAINFODESCRIPTION_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Description of the meta data arrays of MSSpectrum.
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI MetaInfoDescription
    : public MetaInfoInterface
  {
    public:
    	/// Constructor
      MetaInfoDescription();
      /// Copy constructor
      MetaInfoDescription(const MetaInfoDescription& source);
      /// Destructor
      ~MetaInfoDescription();
			
			/// Assignment operator
      MetaInfoDescription& operator= (const MetaInfoDescription& source);

      /// Equality operator
      bool operator== (const MetaInfoDescription& rhs) const;

			/// returns the name of the peak annotations
      const String& getName() const;
      /// sets the name of the peak annotations
      void setName(const String& name);

			/// returns a const reference to the description of the applied processing
      const std::vector<DataProcessing>& getDataProcessing() const;
      /// returns a mutable reference to the description of the applied processing
      std::vector<DataProcessing>& getDataProcessing();
      /// sets the description of the applied processing
      void setDataProcessing(const std::vector<DataProcessing>& data_processing);

    protected:
      String comment_;
      String name_;
      std::vector<DataProcessing> data_processing_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_METAINFODESCRIPTION_H
