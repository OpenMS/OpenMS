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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SOURCEFILE_H
#define OPENMS_METADATA_SOURCEFILE_H

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Description of a file location, used to store the origin of (meta) data.
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI SourceFile
  	: public MetaInfoInterface
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

    	///Native ID type
    	enum NativeIDType
    	{
    		UNKNOWN_NATIVEID,			///< Unknown native ID type
    		THERMO,								///< controller=xsd:nonNegativeInteger scan=xsd:positiveInteger
    		WATERS,								///< function=xsd:positiveInteger process=xsd:nonNegativeInteger scan=xsd:nonNegativeInteger
    		WIFF,									///< sample=xsd:nonNegativeInteger period=xsd:nonNegativeInteger cycle=xsd:nonNegativeInteger experiment=xsd:nonNegativeInteger
    		BRUKER_AGILENT,				///< scan=xsd:nonNegativeInteger
    		BRUKER_BAF,						///< scan=xsd:nonNegativeInteger
    		BRUKER_FID,						///< file=xsd:IDREF
    		BRUKER_U2,            ///< declaration=xsd:nonNegativeInteger collection=xsd:nonNegativeInteger scan=xsd:nonNegativeInteger
    		MULTIPLE_PEAK_LISTS,	///< index=xsd:nonNegativeInteger @n Used for conversion of peak list files with multiple spectra, i.e. MGF, PKL, merged DTA files. Index is the spectrum number in the file, starting from 0.
    		SINGLE_PEAK_LIST,			///< file=xsd:IDREF @n The nativeID must be the same as the source file ID. Used for conversion of peak list files with one spectrum per file, typically folder of PKL or DTAs, each sourceFileRef is different.
    		SCAN_NUMBER,					///< scan=xsd:nonNegativeInteger @n Used for conversion from mzXML, or DTA folder where native scan numbers can be derived.
    		SPECTRUM_IDENTIFIER,	///< spectrum=xsd:nonNegativeInteger @n Used for conversion from mzData. The spectrum id attribute is referenced.
    		AB_SCIEX,             ///< jobRun=xsd:nonNegativeInteger spotLabel=xsd:string spectrum=xsd:nonNegativeInteger
    		AGILENT_MASSHUNTER,   ///< scanId=xsd:nonNegativeInteger
    		SIZE_OF_NATIVEIDTYPE
    	};
			/// Names of native ID types
			static const std::string NamesOfNativeIDType[SIZE_OF_NATIVEIDTYPE];

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
  		NativeIDType getNativeIDType() const;
			/// Sets the native ID type of the spectra
  		void setNativeIDType(NativeIDType type);

    protected:
			String name_of_file_;
			String path_to_file_;
			DoubleReal file_size_;
			String file_type_;
			String checksum_;
			ChecksumType checksum_type_;
    	NativeIDType native_id_type_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_SOURCEFILE_H
