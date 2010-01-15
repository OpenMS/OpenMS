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

#include <OpenMS/METADATA/SourceFile.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	const std::string SourceFile::NamesOfChecksumType[] = {"Unknown","SHA-1","MD5"};

	SourceFile::SourceFile()
		: CVTermList(),
		  name_of_file_(),
		  path_to_file_(),
			file_size_(),
		  file_type_(),
		  checksum_(),
		  checksum_type_(),
    	native_id_type_("")
	{
	  
	}
	
	SourceFile::SourceFile(const SourceFile& source)
		: CVTermList(source),
		  name_of_file_(source.name_of_file_),
		  path_to_file_(source.path_to_file_),
			file_size_(source.file_size_),
		  file_type_(source.file_type_),
		  checksum_(source.checksum_),
		  checksum_type_(source.checksum_type_),
		  native_id_type_(source.native_id_type_)
	{
	}
	
	SourceFile::~SourceFile()
	{
	}
	
	SourceFile& SourceFile::operator = (const SourceFile& source)
	{
	  if (&source == this) return *this;
	  
	  CVTermList::operator=(source);
	  name_of_file_ = source.name_of_file_;
	  path_to_file_ = source.path_to_file_;
	  file_size_ = source.file_size_;
	  file_type_ = source.file_type_;
	  checksum_ = source.checksum_;
	  checksum_type_ = source.checksum_type_;
	  native_id_type_ = source.native_id_type_;
	  
	  return *this;
	}
	
	bool SourceFile::operator== (const SourceFile& rhs) const
	{
		return 
	  	CVTermList::operator==(rhs) &&
	    name_of_file_ == rhs.name_of_file_ &&
	    path_to_file_ == rhs.path_to_file_ &&
		  file_size_ == rhs.file_size_ &&
	    file_type_ == rhs.file_type_ &&
	    checksum_ == rhs.checksum_ &&
	  	checksum_type_ == rhs.checksum_type_ &&
	  	native_id_type_ == rhs.native_id_type_
			;
	}
	
	bool SourceFile::operator!= (const SourceFile& rhs) const
	{
		return !(operator==(rhs));
	}
	
	const String& SourceFile::getNameOfFile() const 
	{
	  return name_of_file_; 
	}
	
	void SourceFile::setNameOfFile(const String& name_of_file)
	{
		name_of_file_ = name_of_file; 
	}
	
	const String& SourceFile::getPathToFile() const 
	{
	  return path_to_file_; 
	}
	
	void SourceFile::setPathToFile(const String& path_to_file)
	{
	  path_to_file_ = path_to_file; 
	}

	Real SourceFile::getFileSize() const
	{
		return file_size_;
	}

	void SourceFile::setFileSize(Real file_size)
	{
		file_size_ = file_size;
	}


	const String& SourceFile::getFileType() const 
	{
	  return file_type_; 
	}
	
	void SourceFile::setFileType(const String& file_type)
	{
	  file_type_ = file_type; 
	}
	
	const String& SourceFile::getChecksum() const
	{
		return checksum_;
	}

	SourceFile::ChecksumType SourceFile::getChecksumType() const
	{
		return checksum_type_;
	}
	
	void SourceFile::setChecksum(const String& checksum, ChecksumType type)
	{
		checksum_ = checksum;
		checksum_type_ = type;
	}

  const String& SourceFile::getNativeIDType() const
	{
		return native_id_type_;
	}
	
  void SourceFile::setNativeIDType(const String& type)
	{
		native_id_type_ = type;
	}

}

