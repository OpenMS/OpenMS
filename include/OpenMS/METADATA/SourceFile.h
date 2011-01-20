// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_METADATA_SOURCEFILE_H
#define OPENMS_METADATA_SOURCEFILE_H

#include <OpenMS/METADATA/CVTermList.h>

namespace OpenMS 
{
	/**
		@brief Description of a file location, used to store the origin of (meta) data.
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI SourceFile
  	: public CVTermList
  {
    public:
    	///Type of the checksum
    	enum ChecksumType
    	{
    		UNKNOWN_CHECKSUM, ///< Unknown checksum type
    		SHA1, 	 ///< Secure Hash Algorithm-1
    		MD5,		 ///< Message-Digest algorithm 5
    		SIZE_OF_CHECKSUMTYPE
    	};
			/// Names of checksum types
			static const std::string NamesOfChecksumType[SIZE_OF_CHECKSUMTYPE];

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
      Real getFileSize() const;
      /// sets the file size in MB
      void setFileSize(Real file_size);

 			/// returns the file type
      const String& getFileType() const;
     	/// sets the file type
      void setFileType(const String& file_type);
      
      /// returns the file's checksum
      const String& getChecksum() const;
      /// sets the file's checksum
      void setChecksum(const String& checksum, ChecksumType type);
      /// returns the checksum type
      ChecksumType getChecksumType() const;

			/// Returns the native ID type of the spectra
  		const String& getNativeIDType() const;
			/// Sets the native ID type of the spectra
  		void setNativeIDType(const String& type);

    protected:
			String name_of_file_;
			String path_to_file_;
			DoubleReal file_size_;
			String file_type_;
			String checksum_;
			ChecksumType checksum_type_;
    	String native_id_type_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_SOURCEFILE_H
