// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/CONCEPT/Exception.h>

#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QStringList>
#include <QtNetwork/QHostInfo>

#include <iostream>

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
	
	String File::absolutePath(const String& file)
	{
		QFileInfo fi(file.c_str());
		return fi.absoluteFilePath();
	}

	String File::basename(const String& file)
	{
		QFileInfo fi(file.c_str());
		return fi.fileName();
	}

	String File::path(const String& file)
	{
		QFileInfo fi(file.c_str());
		return fi.path();
	}

	bool File::readable(const String& file)
	{
		QFileInfo fi(file.c_str());
		return (fi.exists() && fi.isReadable());
	}

	bool File::writable(const String& file)
	{
		QFileInfo fi(file.c_str());
		
		bool tmp(false);
		if (!fi.exists())
		{
			QFile f;
			f.setFileName(file.c_str());
			f.open(QIODevice::WriteOnly);
			tmp = f.isWritable();
			f.close();
		}
		else
		{
			tmp = fi.isWritable();
		}
		
		return tmp;
	}

	String File::find(const String& filename, vector<String> directories)
	{
		String filename_new = filename;
		
		//add env $OPENMS_DATA_PATH
		if (getenv("OPENMS_DATA_PATH") != 0)
		{
			directories.push_back(String(getenv("OPENMS_DATA_PATH")));
		}
		
		//add data dir in OpenMS data path
		directories.push_back(OPENMS_DATA_PATH);
		
		//add path suffix to all specified directories
		String path = File::path(filename);
		if (path!="")
		{
			for (vector<String>::iterator it=directories.begin(); it!=directories.end(); ++it)
			{
				it->ensureLastChar('/');
				*it += path;
			}
			filename_new = File::basename(filename);
		}
		
		//look up file
		for (vector<String>::const_iterator it=directories.begin(); it!=directories.end(); ++it)
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

  bool File::createSparseFile(const String& filename, const Int64& sfilesize = 1)
  { 
		Int64 filesize = sfilesize; 
		if (filesize < 1)
		{
			std::cerr << "File::createSparseFile (__LINE__) warning: mapping of empty files not allowed! Increasing filesize to 1 byte!" << std::endl;
			filesize = 1;
		}
    #ifdef OPENMS_WINDOWSPLATFORM
  
      #ifdef UNICODE
        int len;
        int slength = (int)filename.length() + 1;
        len = MultiByteToWideChar(CP_ACP, 0, filename.c_str(), slength, 0, 0); 
        wchar_t* buf = new wchar_t[len];
        MultiByteToWideChar(CP_ACP, 0, filename.c_str(), slength, buf, len);
        std::wstring stemp(buf);
        delete[] buf;
        LPCWSTR result = stemp.c_str();
      #else
        LPCTSTR result = filename.c_str();
      #endif
  
      HANDLE hFile = CreateFile( result , GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
      if (hFile == INVALID_HANDLE_VALUE)
      {
        return false;
      }
      DWORD dwTemp;
      if (DeviceIoControl(hFile, FSCTL_SET_SPARSE, NULL, 0, NULL, 0, &dwTemp, NULL) == 0)
      {
        CloseHandle(hFile);
        return false;
      }
      LARGE_INTEGER fs;
      fs.QuadPart = filesize;
      if (SetFilePointerEx(hFile, fs, NULL, FILE_BEGIN) == 0)
      {
        CloseHandle(hFile);
        return false;
      }
      if (SetEndOfFile(hFile) == 0)
      {
        CloseHandle(hFile);
        return false;
      }
      CloseHandle(hFile);
    
    #else
    
      int fd = open64(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
      if (fd == -1) {
				std::cerr << "failed while opening " << filename << "\n";
        return false;
      }
      
      /* Stretch the file size
      */
			//std::cerr << "!!Seeking to "<< filesize << " in " << filename << "\n";
      int result = lseek64(fd, filesize, SEEK_SET);
      if (result == -1) {
				std::cerr << "failed while seeking to "<< filesize << " in " << filename << "\n";
        close(fd);
        return false;
      }
      
      /* Something needs to be written at the end of the file to
      * have the file actually have the new size.
      * Just writing an empty string at the current file position will do.
      */
      result = write(fd, "", 1);
      close(fd);
      if (result != 1) {
				std::cerr << "failed while writing to " << filename << "\n";
        return false;
      }
    #endif
    return true;
  }

	#ifdef OPENMS_WINDOWSPLATFORM
  bool File::extendSparseFile(const HANDLE& /*hFile*/, const Int64& /*filesize*/)
	#else
  bool File::extendSparseFile(const int& hFile, const Int64& filesize)
  #endif
  { 
		//TODO see http://www.boost.org/libs/filesystem/doc/tr2_proposal.html#Text 
		//     for template <class Path> space_info space(const Path& p);
		
		#ifdef OPENMS_WINDOWSPLATFORM
		// for Windows, this is a null-operation, as the OS will extend memory mapped files automatically
		// as long as the call to 'CreateFileMapping' allows for a large enough mapping
		return true;
		
		#else
      
		/* Stretch the file size
		*/
		int result = lseek64(hFile, filesize, SEEK_SET);
		if (result == -1) {
			std::cerr << "failed while seeking to "<< filesize << " in extendSparseFile" << "\n";
			return false;
		}
		
		/* Something needs to be written at the end of the file to
		* have the file actually have the new size.
		* Just writing an empty string at the current file position will do.
		*/
		result = write(hFile, "", 1);
		if (result != 1) {
			std::cerr << "failed while writing in extendSparseFile\n";
			return false;
		}
		return true;
		#endif
  }	
	
  #ifdef OPENMS_WINDOWSPLATFORM
  HANDLE File::getSwapFileHandle(const String& filename, const Int64& filesize, const bool& create)
  {
    if (create && (!File::exists(filename)))
    {
      if (!createSparseFile( filename, filesize ))
      {
        throw Exception::UnableToCreateFile( __FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToCreateFile in getSwapFileHandle");
      }
    }
    // create mapping object (needed for windows-mmap call)
    #ifdef UNICODE
      int len;
      int slength = (int)filename.length() + 1;
      len = MultiByteToWideChar(CP_ACP, 0, filename.c_str(), slength, 0, 0); 
      wchar_t* buf = new wchar_t[len];
      MultiByteToWideChar(CP_ACP, 0, filename.c_str(), slength, buf, len);
      std::wstring stemp(buf);
      delete[] buf;
      LPCWSTR result = stemp.c_str();
    #else
      LPCTSTR result = filename.c_str();
    #endif
    HANDLE myFile = CreateFile( result , 
                                FILE_WRITE_DATA | FILE_READ_DATA, 
                                FILE_SHARE_DELETE|FILE_SHARE_READ|FILE_SHARE_WRITE, 
                                NULL, 
                                OPEN_EXISTING,
                                FILE_ATTRIBUTE_TEMPORARY, //alternative: FILE_ATTRIBUTE_NORMAL
                                NULL);
                                
		if (myFile == INVALID_HANDLE_VALUE)
		{
			throw Exception::FileNotFound( __FILE__, __LINE__, __PRETTY_FUNCTION__, filename.c_str());
		}
                                
    return myFile;
  }
  #else
  int File::getSwapFileHandle(const String& filename, const Int64& filesize, const bool& create)
  {
    if (create && (!File::exists(filename)))
    {
      if (!createSparseFile( filename, filesize ))
      {
        throw Exception::UnableToCreateFile( __FILE__, __LINE__, __PRETTY_FUNCTION__, filename.c_str());
      }
    }
    int mmapHandle_ = open64(filename.c_str(), O_RDWR, (mode_t)0600);
    if (mmapHandle_ == -1)
    {
      throw Exception::FileNotFound( __FILE__, __LINE__, __PRETTY_FUNCTION__, filename.c_str());
    }
    return mmapHandle_;
  }        
  #endif  

  #ifdef OPENMS_WINDOWSPLATFORM
  void File::closeSwapFileHandle(const HANDLE & f_handle)
  {
    CloseHandle(f_handle);
  }
  #else
  void File::closeSwapFileHandle(const int & f_handle)
  {
    close(f_handle);
  }
  #endif    
    
  



  
} // namespace OpenMS
