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

#include <OpenMS/METADATA/SourceFile.h>

#include<iostream>

using namespace std;

namespace OpenMS
{
	SourceFile::SourceFile():
	  name_of_file_(),
	  path_to_file_(),
		file_size_(),
	  file_type_(),
	  sha1_()
	{
	  
	}
	
	SourceFile::SourceFile(const SourceFile& source):
	  name_of_file_(source.name_of_file_),
	  path_to_file_(source.path_to_file_),
		file_size_(source.file_size_),
	  file_type_(source.file_type_),
	  sha1_(source.sha1_)
	{
	  
	}
	
	SourceFile::~SourceFile()
	{
	  
	}
	
	SourceFile& SourceFile::operator = (const SourceFile& source)
	{
	  if (&source == this) return *this;
	  
	  name_of_file_ = source.name_of_file_;
	  path_to_file_ = source.path_to_file_;
	  file_size_ = source.file_size_;
	  file_type_ = source.file_type_;
	  sha1_ = source.sha1_;
	  
	  return *this;
	}
	
	bool SourceFile::operator== (const SourceFile& rhs) const
	{
		return 
	    name_of_file_ == rhs.name_of_file_ &&
	    path_to_file_ == rhs.path_to_file_ &&
		  file_size_ == rhs.file_size_ &&
	    file_type_ == rhs.file_type_ &&
	    sha1_ == rhs.sha1_
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

	const float& SourceFile::getFileSize() const
	{
		return file_size_;
	}

	void SourceFile::setFileSize(const float& file_size)
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
	
	const String& SourceFile::getSha1() const
	{
		return sha1_;
	}
	
	void SourceFile::setSha1(const String& sha1)
	{
		sha1_ = sha1;
	}
	
	bool SourceFile::isFileEmpty() const
	{
    return
    (
    	name_of_file_.size() == 0 &&
    	path_to_file_.size() == 0 &&
    	file_size_ == 0 &&
    	sha1_ == "" &&
    	file_type_ == ""
    );
	}
}

