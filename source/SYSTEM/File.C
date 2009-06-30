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

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QStringList>
#include <QtNetwork/QHostInfo>

#include <iostream>
#include <cstdio>

#ifdef OPENMS_WINDOWSPLATFORM
#  include <Winioctl.h> // for DeviceIoControl and constants e.g. FSCTL_SET_SPARSE
#else
#  include <fcntl.h> // for O_RDWR etc
#  include <unistd.h>
#endif

// Mac OS X does not provide lseek64 and open64, so we need to replace them with their normal counterparts
#if defined __APPLE__ & defined __MACH__
#define lseek64 lseek
#define open64  open 
#endif


using namespace std;

namespace OpenMS 
{
	
	bool File::exists(const String& file)
	{
		QFileInfo fi(file.toQString());
		return fi.exists();
	}

	bool File::empty(const String& file)
	{
		QFileInfo fi(file.toQString());
		return (!fi.exists() || fi.size()==0);
	}

	bool File::remove(const String& file)
	{
		if (!exists(file)) return true;
		
	  if( std::remove(file.c_str()) != 0 ) return false;
	  return true;
	}
	
	String File::absolutePath(const String& file)
	{
		QFileInfo fi(file.toQString());
		return fi.absoluteFilePath();
	}

	String File::basename(const String& file)
	{
		QFileInfo fi(file.toQString());
		return fi.fileName();
	}

	String File::path(const String& file)
	{
		QFileInfo fi(file.toQString());
		return fi.path();
	}

	bool File::readable(const String& file)
	{
		QFileInfo fi(file.toQString());
		return (fi.exists() && fi.isReadable());
	}

	bool File::writable(const String& file)
	{
		QFileInfo fi(file.toQString());
		
		bool tmp = false ;
		if (fi.exists())
		{
			tmp = fi.isWritable();
		}
		else
		{
			QFile f;
			f.setFileName(file.toQString());
			f.open(QIODevice::WriteOnly);
			tmp = f.isWritable();
			f.remove();
		}
		
		return tmp;
	}

	String File::find(const String& filename, StringList directories)
	{
		String filename_new = filename;
		
		//add data dir in OpenMS data path
		directories.push_back(getOpenMSDataPath());
		
		//add path suffix to all specified directories
		String path = File::path(filename);
		if (path!="")
		{
			for (StringList::iterator it=directories.begin(); it!=directories.end(); ++it)
			{
				it->ensureLastChar('/');
				*it += path;
			}
			filename_new = File::basename(filename);
		}
		
		//look up file
		for (StringList::const_iterator it=directories.begin(); it!=directories.end(); ++it)
		{
			String loc = *it;
			loc.ensureLastChar('/');
			loc = loc + filename_new;
			
			if (exists(loc))
			{
				return loc;
			}
		}
		
		//if the file was not found, throw an exception
		throw Exception::FileNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,filename);
		
		//this is never reached, but needs to be there to avoid compiler warnings
		return "";
	}
	
	bool File::fileList(const String& dir, const String& file_pattern, StringList& output)
	{
		QDir d(dir.toQString(), file_pattern.toQString(), QDir::Name, QDir::Files);
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
			output[i++] = (*it);
		}
		
		return true;
	}

	String File::getUniqueName()
	{
		DateTime now = DateTime::now();
		String pid;
		#ifdef OPENMS_WINDOWSPLATFORM
			pid = (String)GetCurrentProcessId();
		#else
			pid = (String)getpid();	
		#endif		
		static int number = 0;
		return now.getDate() + "_" + now.getTime().remove(':') + "_" + String(QHostInfo::localHostName()) + "_" + pid + "_" + (++number);
	}
  
  String File::getOpenMSDataPath()
  {
		if (getenv("OPENMS_DATA_PATH") != 0)
		{
			return getenv("OPENMS_DATA_PATH");
		}
		
		return OPENMS_DATA_PATH;
  }

	String File::removeExtension(const OpenMS::String& file) 
	{
		if (!file.has('.')) return file;
		
		SignedSize ext_length = file.suffix('.').size() + 1;
		return file.substr(0, - ext_length);
	}

	bool File::isDirectory(const String& path)
	{
		QFileInfo fi(path.toQString());
		return fi.isDir();
	}

} // namespace OpenMS
