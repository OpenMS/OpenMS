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

#include <qfileinfo.h>

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

using namespace std;

namespace OpenMS 
{

	bool File::exists(const string& file)
	{
		QFileInfo fi(file);
		return fi.exists();
	}

	bool File::empty(const string& file)
	{
		QFileInfo fi(file);
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
		QFileInfo fi(file);
		file = fi.absFilePath().ascii();
	}

	bool File::readable(const string& file)
	{
		QFileInfo fi(file);
		return (fi.exists() && fi.isReadable());
	}

	bool File::writable(const string& file)
	{
		QFile f;
		f.setName(file);
		f.open(IO_WriteOnly);
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

} // namespace OpenMS
