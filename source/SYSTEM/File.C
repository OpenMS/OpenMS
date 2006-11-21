// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <qfileinfo.h>

#include <fstream>

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

} // namespace OpenMS
