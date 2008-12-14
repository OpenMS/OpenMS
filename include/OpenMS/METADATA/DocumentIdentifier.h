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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_DOCUMENTIDENTIFIER_H
#define OPENMS_METADATA_DOCUMENTIDENTIFIER_H

// OpenMS
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /** 
    @brief Manage source document information.
    
    This class stored information about the source document.
    Primarily this is the document id e.g. a LSID.
    
    For source files additional information can be stored:
    - file name
    - file type
    
    @todo Add source file name and type (FileHandler::Type). Fill these members when loading a map from file (featureXML, consensusXML, mzData, mzXML, mzML, ANDI/MS, DTA, DTA2D, MGF). Use this informtation e.g. when converting to mzML where it is required (Hiwi)
    
    @ingroup Metadata
  */
  
  class OPENMS_DLLAPI DocumentIdentifier
  {
    public:
      /** @name Constructors and Destructors
      */
      //@{
      /// default constructor
      DocumentIdentifier();

      /// Copy constructor
      DocumentIdentifier(const DocumentIdentifier& source);

      /// Assignment operator
      DocumentIdentifier& operator=(const DocumentIdentifier& source);

			/// Equality operator
			bool operator== (const DocumentIdentifier& rhs) const;
			
      /// destructor
      virtual ~DocumentIdentifier();
      //@}    

      /** @name Acessors
       */
      //@{

			/// retrieve computed zero-charge feature map
      void setIdentifier(const String& id);

      /// retrieve computed zero-charge feature map
      const String& getIdentifier() const;

			/// exchange content with @p from
			void swap(DocumentIdentifier& from);
			
      //@}
      
    protected:
      /// the ID (e.g. LSID)
      String id_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_DOCUMENTIDENTIFIER_H
