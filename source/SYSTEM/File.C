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

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QStringList>

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

using namespace std;

namespace OpenMS 
{

	bool File::exists(const string& file)
	{
		QFileInfo fi(file.c_str());
		return fi.exists();
	}

	bool File::empty(const string& file)
	{
		QFileInfo fi(file.c_str());
		return (!fi.exists() || fi.size()==0);
	}

	bool File::remove(const string& file)
	{
		if (!exists(file)) return true;
		
	  if( std::remove(file.c_str()) != 0 ) return false;
	  return true;
	}
	
	void File::absolutePath(string& file)
	{
		QFileInfo fi(file.c_str());
		file = fi.absoluteFilePath().toAscii().data();
	}

	String File::basename(const string& file)
	{
		QFileInfo fi(file.c_str());
		return fi.fileName().toAscii().data();
	}

	String File::path(const string& file)
	{
		QFileInfo fi(file.c_str());
		return fi.filePath().toAscii().data();
	}

	bool File::readable(const string& file)
	{
		QFileInfo fi(file.c_str());
		return (fi.exists() && fi.isReadable());
	}

	bool File::writable(const string& file)
	{
		QFile f;
		f.setFileName(file.c_str());
		f.open(QIODevice::WriteOnly);
		bool tmp = f.isWritable();
		f.close();
		
		return tmp;
	}

	String File::find(const String& filename, vector<String> directories)
	{
		//add env $OPENMS_PATH
		if (getenv("OPENMS_PATH") != 0)
		{
			directories.push_back(String(getenv("OPENMS_PATH")) + "/data/");
		}
		
		//add data dir in OpenMS built path
		directories.push_back(OPENMS_PATH"/data/");
		
		//look up file
		for (vector<String>::const_iterator it=directories.begin(); it!=directories.end(); ++it)
		{
			String loc = *it;
			loc.ensureLastChar('/');
			loc = loc + filename;
			if (exists(loc))
			{
				return loc;
			}
		}
		
		return "";
	}
	
	bool File::fileList(const std::string& dir, const std::string& file_pattern, vector<String>& output)
	{
		QDir d(dir.c_str(), file_pattern.c_str(), QDir::Name, QDir::Files);
		QStringList list = d.entryList();

		//clear and check if empty
		output.clear();
		if (list.size()==0)
		{
			return false;
		}
		
		//resize output
		output.resize(list.size());
		
		//fill output
		UnsignedInt i = 0;
		for ( QStringList::const_iterator it = list.constBegin(); it != list.constEnd(); ++it )
		{
			output[i++] = (*it).toAscii().data();
		}
		
		return true;
	}

} // namespace OpenMS
