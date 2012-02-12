// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_DOCUMENTIDENTIFIER_H
#define OPENMS_METADATA_DOCUMENTIDENTIFIER_H

// OpenMS
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FileTypes.h>

namespace OpenMS
{
  /**
     @brief Manage source document information.

     This class stored information about the source document.
     Primarily this is the document id e.g. a LSID.

     For source files additional information can be stored:
     - file name
     - file type

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

			/// set document identifier (e.g. an LSID)
      void setIdentifier(const String& id);

      /// retrieve document identifier (e.g. an LSID)
      const String& getIdentifier() const;

			/// exchange content with @p from
			void swap(DocumentIdentifier& from);


      /// set the file_name_ according to absolute path of the file loaded from preferrably done whilst loading
      void setLoadedFilePath(const String& file_name);

      /// get the file_name_ which is the absolute path to the file loaded from
      const String& getLoadedFilePath() const;

      /// set the file_type according to the type of the file loaded from (see FileHandler::Type) preferrably done whilst loading
      void setLoadedFileType(const String& file_name);

      /// get the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded from
      const FileTypes::Type& getLoadedFileType() const;

      //@}

    protected:
      /// the ID (e.g. LSID)
      String id_;
      /// the path to the loaded file
      String file_path_;
      /// the type of the loaded file
      FileTypes::Type file_type_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_DOCUMENTIDENTIFIER_H
