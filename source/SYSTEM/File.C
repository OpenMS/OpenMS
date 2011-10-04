// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QStringList>
#include <QtNetwork/QHostInfo>

#include <iostream>
#include <cstdio>

#ifdef OPENMS_WINDOWSPLATFORM
#  define NOMINMAX
#  include <Windows.h>  // for GetCurrentProcessId() && GetModuleFileName()
#else
#  include <unistd.h>   // for 'getpid()'
#endif

using namespace std;

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

namespace OpenMS 
{
	
  String File::getExecutablePath()
  {
    // see http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe/1024937#1024937 for more OS' (if needed)
    static String spath = "";
    static bool path_checked = false;

    if (path_checked) return spath; // short route. Only inquire the path once. The result will be the same every time.

    char path[1024];

#ifdef OPENMS_WINDOWSPLATFORM
    int size = sizeof(path);
    if( GetModuleFileName( NULL, path, size ) )
#elif  defined(__APPLE__)
    uint size = sizeof(path);
    if (_NSGetExecutablePath(path, &size) == 0)
#else // LINUX
    int size = sizeof(path);
    int ch = readlink("/proc/self/exe", path, size);
    if (ch != -1)
#endif
    { 
      spath = File::path(String(path));
      if (File::exists(spath)) // check if directory exists
      {
        // ensure path ends with a "/", such that we can just write path + "ToolX", and to not worry about if its empty or a path.
        spath.ensureLastChar('/');
      }
      else
      {
        std::cerr << "Path extracted from Executable Path does not exist! Returning empty string!\n";
        spath = "";
      }
    }
    else
    {
      std::cerr << "Cannot get Executable Path! Not using a path prefix!\n";
    }

    path_checked = true; // enable short route for next run
    return spath;
  }

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
	
	bool File::removeDirRecursively(const String& dir_name)
	{
		bool fail = false;
		QString path = dir_name.toQString();
    QDir dir(path);
		QStringList files = dir.entryList(QDir::Files | QDir::NoDotAndDotDot);
		foreach (const QString& file_name, files)
		{
			if (!dir.remove(file_name))
			{
				LOG_WARN << "Could not remove file " << String(file_name) << "!" << std::endl;
				fail = true;
			}
		}
		QStringList contained_dirs = dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
		foreach (const QString& contained_dir, contained_dirs)
		{
			if (!removeDirRecursively(path + QDir::separator() + contained_dir))
			{
				fail = true;
			}
		}
		
		QDir parent_dir(path);
		if (parent_dir.cdUp())
		{
			if (!parent_dir.rmdir(path))
			{
				std::cerr << "Could not remove directory " << String(dir.dirName()) << "!" << std::endl;
				fail = true;
			}
		}
		
		return !fail;
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
        return String(QDir::cleanPath(loc.toQString()));
			}
		}
		
		//if the file was not found, throw an exception
		throw Exception::FileNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,filename);
		
		//this is never reached, but needs to be there to avoid compiler warnings
		return "";
	}

	bool File::fileList(const String& dir, const String& file_pattern, StringList& output, bool full_path)
	{
		QDir d(dir.toQString(), file_pattern.toQString(), QDir::Name, QDir::Files);
		QFileInfoList list = d.entryInfoList();

		//clear and check if empty
		output.clear();
		if (list.empty())
		{
			return false;
		}
		
		//resize output
		output.resize(list.size());
		
		//fill output
		UInt i = 0;
		for (QFileInfoList::const_iterator it = list.constBegin(); it != list.constEnd(); ++it)
		{
			output[i++] = full_path ? it->filePath() : it->fileName();
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
		String path;
		if (getenv("OPENMS_DATA_PATH") != 0)
		{
			path = getenv("OPENMS_DATA_PATH");
		}
		else
		{
			path = OPENMS_DATA_PATH;
		}
		
		if (!exists(path))
		{ // now we're in big trouble as './share' is not were its supposed to be...
			std::cerr << "OpenMS FATAL ERROR!\nExpected shared data to be at '" << path << "'! OpenMS cannot function without it! Exiting now. To resolve this, set the environment variable 'OPENMS_DATA_PATH' to the OpenMS share directory.\n";
			exit(1);
		}
		
		return path;
  }

	String File::removeExtension(const OpenMS::String& file) 
	{
		if (!file.has('.')) return file;
		
		SignedSize ext_length = file.suffix('.').size() + 1;
		return file.chop(ext_length);
	}

	bool File::isDirectory(const String& path)
	{
		QFileInfo fi(path.toQString());
		return fi.isDir();
	}

  String File::getTempDirectory()
  {
    Param p = getSystemParameters();
    if (p.exists("temp_dir") && String(p.getValue("temp_dir")).trim() != "")
    {
       return p.getValue("temp_dir");
    }
    return String(QDir::tempPath());
  }

  /// The current OpenMS user data path (for result files)
  String File::getUserDirectory()
  {
    Param p = getSystemParameters();
    String dir;
    if (p.exists("home_dir") && String(p.getValue("home_dir")).trim() != "")
    {
      dir = p.getValue("home_dir");
    }
    else
    {
      dir = String(QDir::homePath());
    }
    dir.ensureLastChar('/');
    return dir;
  }

  String File::findDatabase(const String& db_name)
  {
    Param sys_p = getSystemParameters();
    String full_db_name;
    try
    {
      full_db_name = find(db_name, sys_p.getValue("id_db_dir"));
      LOG_INFO << "Augmenting database name '" << db_name << "' with path given in 'OpenMS.ini:id_db_dir'. Full name is now: '" << full_db_name << "'\n";
    }
    catch (Exception::BaseException& e)
    {
      LOG_ERROR << "Input database '" + db_name + "' not found. Make sure it exists (and check 'OpenMS.ini:id_db_dir' if you used relative paths. Aborting!";
      throw e;
    }

    return full_db_name;
  }

  Param File::getSystemParameters()
  {
    String filename = String(QDir::homePath()) + "/.OpenMS/OpenMS.ini";
    Param p;
    if (!File::readable(filename))
    { // create file
      p = getSystemParameterDefaults_();
      QDir qd;
      if (!qd.exists(String(String(QDir::homePath()) + "/.OpenMS/").toQString()))
      {
        qd.mkpath(String(String(QDir::homePath()) + "/.OpenMS/").toQString());
      }
      p.store(filename);
    }
    else
    {
      p.load(filename);
      // check version
      if (!p.exists("version") || (p.getValue("version") != VersionInfo::getVersion()))
      {
			  if (!p.exists("version"))
        {
          LOG_WARN << "Broken file '"<< filename << "' discovered. The 'version' tag is missing.\n";
        }
        else // old version
        {
          LOG_WARN << "File '"<< filename << "' is deprecated.\n";
        }
        LOG_WARN << "Updating missing/wrong entries in '"<< filename << "' with defaults!\n";
        Param p_new = getSystemParameterDefaults_();
        p.setValue("version", VersionInfo::getVersion()); // update old version, such that p_new:version does not get overwritten during update()
        p_new.update(p);
        p_new.store(filename);
      }
    }
    return p;
  }

  Param File::getSystemParameterDefaults_()
  {
    Param p;
    p.setValue("version", VersionInfo::getVersion());
    p.setValue("home_dir", ""); // only active when user enters something in this value
    p.setValue("temp_dir", ""); // only active when user enters something in this value
    p.setValue("id_db_dir", StringList::create(""),
                            String("Default directory for FASTA and psq files used as databased for id engines. ") + \
                                   "This allows you to specify just the filename of the DB in the " + \
                                   "respective TOPP tool, and the database will be searched in the directories specified here " + \
                                   "");// only active when user enters something in this value
    p.setValue("threads", 1);
    // TODO: maybe we add -log, -debug.... or....

    return p;
  }


} // namespace OpenMS
