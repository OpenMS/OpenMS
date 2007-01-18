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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SOURCEFILE_H
#define OPENMS_METADATA_SOURCEFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS 
{
	/**
		@brief Description of a file location, used to store the origin of (meta) data.
		
		
		
		@ingroup Metadata
	*/
  class SourceFile
  {
    public:
    	/// Constructor
      SourceFile();
      /// Copy constructor
      SourceFile(const SourceFile& source);
      /// Destructor
      ~SourceFile();
      /// Assignment operator
      SourceFile& operator= (const SourceFile& source);

      /// Equality operator
      bool operator== (const SourceFile& rhs) const;
      /// Equality operator
      bool operator!= (const SourceFile& rhs) const;
			
			/// returns the file name
      const String& getNameOfFile() const;
      /// sets the file name
      void setNameOfFile(const String& name_of_file);

			/// returns the file path
      const String& getPathToFile() const;
      /// sets the file path
      void setPathToFile(const String& path_path_to_file);
			
      /// returns the file size in MB
      const float& getFileSize() const;
      /// sets the file size in MB
      void setFileSize(const float& file_size);

 			/// returns the file type
      const String& getFileType() const;
     	/// sets the file type
      void setFileType(const String& file_type);
      
      /// returns the source file's SHA1 hash value
      const String& getSha1() const;
      /// sets the source file's SHA1 hash value
      void setSha1(const String& sha1);

	    /// returns if the SourceFile is empty
	    bool isFileEmpty() const;

    protected:
			String name_of_file_;
			String path_to_file_;
			float file_size_;
			String file_type_;
			String sha1_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_SOURCEFILE_H
