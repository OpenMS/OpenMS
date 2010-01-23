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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/FORMAT/FileHandler.h>


namespace OpenMS
{

	DocumentIdentifier::DocumentIdentifier()
		: id_(),
			file_path_(),
			file_type_(FileTypes::UNKNOWN)
	{
	}

	DocumentIdentifier::DocumentIdentifier(const DocumentIdentifier& source)
		: id_(source.id_),
	    file_path_(source.file_path_),
			file_type_(source.file_type_)
	{
	}

	DocumentIdentifier& DocumentIdentifier::operator=(const DocumentIdentifier& source)
	{
		if (&source==this)
		{
			return *this;
		}

		id_ = source.id_;
		file_type_ = source.file_type_;
		file_path_ = source.file_path_;

		return *this;
	}

	DocumentIdentifier::~DocumentIdentifier()
	{
	}

	void DocumentIdentifier::setIdentifier(const String& id)
	{
		id_ = id;
	}

	const String& DocumentIdentifier::getIdentifier() const
	{
		return id_;
	}

	void DocumentIdentifier::setLoadedFilePath(const String& file_name)
	{
		file_path_ = File::absolutePath(file_name);
	}

	const String& DocumentIdentifier::getLoadedFilePath() const
	{
		return file_path_;
	}

	void DocumentIdentifier::setLoadedFileType(const String& file_name)
	{
		file_type_ = FileHandler::getTypeByContent(file_name);
	}

	const FileTypes::Type& DocumentIdentifier::getLoadedFileType() const
	{
		return file_type_;
	}

	void DocumentIdentifier::swap(DocumentIdentifier& from)
	{
		std::swap(id_, from.id_);
		std::swap(file_path_, from.file_path_);
		std::swap(file_type_, from.file_type_);
	}

	bool DocumentIdentifier::operator== (const DocumentIdentifier& rhs) const
  {
  	return  id_ == rhs.id_;
	}


}

