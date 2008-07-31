// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_METAINFODESCRIPTION_H
#define OPENMS_METADATA_METAINFODESCRIPTION_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Description of the meta data arrays of DSpectrum.
		
		@ingroup Metadata
	*/
  class MetaInfoDescription
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
			
			/// returns the free-text comment
      const String& getComment() const;
      /// sets the free-text comment
      void setComment(const String& comment);
			
			/// returns a const reference to the source file the information
      const SourceFile& getSourceFile() const;
      /// returns a mutable reference to the source file the information
      SourceFile& getSourceFile();
      /// sets the source file the information
      void setSourceFile(const SourceFile& source_file);

			/// returns the name of the peak annotations
      inline const String& getName() const
      {
      	return name_;
      }
      /// sets the name of the peak annotations
      void setName(const String& name);
      
    protected:
      String comment_;
      String name_;
      SourceFile source_file_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_METAINFODESCRIPTION_H
