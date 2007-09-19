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
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QStringList>
#include <QtNetwork/QHostInfo>

#ifdef OPENMS_OS_MINGW32
#include <windows.h>
#endif

using namespace std;

namespace OpenMS 
{
	
	bool File::exists(const String& file)
	{
		QFileInfo fi(file.c_str());
		return fi.exists();
	}

	bool File::empty(const String& file)
	{
		QFileInfo fi(file.c_str());
		return (!fi.exists() || fi.size()==0);
	}

	bool File::remove(const String& file)
	{
		if (!exists(file)) return true;
		
	  if( std::remove(file.c_str()) != 0 ) return false;
	  return true;
	}
	
	void File::absolutePath(String& file)
	{
		QFileInfo fi(file.c_str());
		file = fi.absoluteFilePath().toAscii().data();
	}

	String File::basename(const String& file)
	{
		QFileInfo fi(file.c_str());
		return fi.fileName().toAscii().data();
	}

	String File::path(const String& file)
	{
		QFileInfo fi(file.c_str());
		return fi.path().toAscii().data();
	}

	bool File::readable(const String& file)
	{
		QFileInfo fi(file.c_str());
		return (fi.exists() && fi.isReadable());
	}

	bool File::writable(const String& file)
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
	
	bool File::fileList(const String& dir, const String& file_pattern, vector<String>& output)
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
		UInt i = 0;
		for ( QStringList::const_iterator it = list.constBegin(); it != list.constEnd(); ++it )
		{
			output[i++] = (*it).toAscii().data();
		}
		
		return true;
	}

	String File::getUniqueName()
	{
		DateTime now;
		String date_str, time_str, pid;
		now.now();
		now.getDate(date_str);
		now.getTime(time_str);
		#ifdef OPENMS_OS_MINGW32
			pid = (String)GetCurrentProcessId();
		#else
			pid = (String)getpid();	
		#endif		
		time_str.remove(':'); // remove ':', because of Windoze 
		return date_str + "_" + time_str + "_" + String(QHostInfo::localHostName()) + "_" + pid;
	}

} // namespace OpenMS
