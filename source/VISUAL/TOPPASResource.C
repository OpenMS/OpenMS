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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <iostream>
#include <OpenMS/VISUAL/TOPPASResource.h>

namespace OpenMS
{
	QStringList TOPPASResource::supported_schemes = (QStringList() << "file");
	
	TOPPASResource::TOPPASResource(const QString& file)
		: QObject(),
			url_(),
			file_name_("")
	{
		fromLocalFile(file);
	}
	
	TOPPASResource::TOPPASResource(const QUrl& url)
		: QObject(),
			url_(),
			file_name_("")
	{
		QString scheme = url.scheme().toLower();
		if (!supported_schemes.contains(scheme))
		{
			std::cerr << "URL scheme not supported!" << std::endl;
		}
		else
		{
			url_ = url;
			
			if (scheme == "file")
			{
				file_name_ = url.toLocalFile();
			}
		}
	}
	
	TOPPASResource::TOPPASResource(const TOPPASResource& rhs)
		: QObject(),
			url_(rhs.url_),
			file_name_(rhs.file_name_)
	{
	}
	
	TOPPASResource::~TOPPASResource()
	{
	}
	
	TOPPASResource& TOPPASResource::operator= (const TOPPASResource& rhs)
	{
		url_ = rhs.url_;
		file_name_ = rhs.file_name_;
		
		return *this;
	}
	
	void TOPPASResource::writeToFile(const QString& file_name)
	{
		// TODO retrieve data and write it to file_name
		
		file_name_ = file_name;
	}
	
	const QString& TOPPASResource::getLocalFile() const
	{
		return file_name_;
	}
	
	const QUrl& TOPPASResource::getURL() const
	{
		return url_;
	}
	
	void TOPPASResource::fromLocalFile(const QString& file)
	{
		url_ = QUrl::fromLocalFile(file);
		file_name_ = file;
	}
	
} //namespace


