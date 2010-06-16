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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/DocumentIDTagger.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDir>

#ifdef _MSC_VER // disable some boost warnings that distract from ours
#	pragma warning( push ) // save warning state
#	pragma warning( disable : 4018 )
#endif
#include <boost/interprocess/sync/file_lock.hpp>
#ifdef _MSC_VER
#	pragma warning( pop )  // restore old warning state
#endif


#include <algorithm>
#include <fstream>
#include <iostream>
#include <exception>
#include <cstdio>
#include <ctime>

using namespace std;

namespace OpenMS 
{
  DocumentIDTagger::DocumentIDTagger(String toolname): 
    toolname_(toolname),
		pool_file_()
  {
		pool_file_ = File::getOpenMSDataPath() + ("/IDPool/IDPool.txt");
  }
  
  DocumentIDTagger::DocumentIDTagger(const DocumentIDTagger& source):
    toolname_(source.toolname_),
	  pool_file_(source.pool_file_)
  {
  }
   
  DocumentIDTagger::~DocumentIDTagger()
  {
  }
  
  DocumentIDTagger& DocumentIDTagger::operator = (const DocumentIDTagger& source)
  {
    if (source == *this) return *this;
		toolname_ = source.toolname_;
		pool_file_ = source.pool_file_;
    return *this;
  }

  bool DocumentIDTagger::operator == (const DocumentIDTagger& rhs) const
  {
    return ( (toolname_ == rhs.toolname_)
					&& (pool_file_ == rhs.pool_file_));
  }

  bool DocumentIDTagger::operator != (const DocumentIDTagger& rhs) const
  {
    return !(operator == (rhs));
  }

	String DocumentIDTagger::getPoolFile() const
	{
		return pool_file_;
	}

	void DocumentIDTagger::setPoolFile(const String& file)
	{
		pool_file_ = file;
	}

	bool DocumentIDTagger::getID_(String& id, Int& free, bool idcount_only) const
	{
		free = 0;

		String IDPool_file = getPoolFile();
		String IDPool_file_tmp = String(IDPool_file) + String(".tmp");
		// create PoolFile if non-existant
		if (!File::exists(IDPool_file)) 
		{
			ofstream out(IDPool_file.c_str());
			out.close();
		}

		// open input file
		ifstream in(IDPool_file.c_str());
		if( !in.is_open())
		{
			std::cerr << "IDTagger::getID_() " << IDPool_file << " file failed to open.\n";
			return false;
		}
	
		// CREATE lock file!
		// Hint: we create an independent lock file, because to remove one ID the pool is copied to another file 
		//       which overwrites the original pool file. Additionally we want to atomically write to a log file.
		//       So locking the pool file itself is a really bad idea!
		String tmp_lock_file = String(IDPool_file + String(".lck"));
		if (!File::exists(tmp_lock_file)) 
		{
			ofstream out(tmp_lock_file.c_str());
			out.close();
		}
		boost::interprocess::file_lock flock(tmp_lock_file.c_str());

		// this might throw an exception!
		try {	flock.lock();} catch (exception /*e*/) {return false;}
		// we have the lock!

		// now open temp output file
		ofstream out;
		if (!idcount_only) 
		{
			out.open(IDPool_file_tmp.c_str(), ios::out | ios::trunc);
			if( !out.is_open())
			{
				std::cerr << "IDTagger::getID_() " << IDPool_file_tmp << " file failed to open for writing.\n";
				flock.unlock();
				in.close();
				return false;
			}
		}

		// copy the ID pool, but excise the first ID
		string line;
		while (! in.eof() )
    {
      getline (in,line);
			if (line.length()==0) continue;
			++free;
			if (free==1) id = line;// pull out first ID
			if (!idcount_only)
			{
				if (free!=1) // delete first line
				{
					out << line << "\n";
				}
			}
		}
		in.close();
		if (!idcount_only) 
		{
			out.close();    
			// delete the original file
			remove(IDPool_file.c_str());
			// rename old to new
			rename(IDPool_file_tmp.c_str(),IDPool_file.c_str());

			// write ID to log file
			String logging_file = String(IDPool_file + String(".log"));
			ofstream outfile;
			outfile.open (logging_file.c_str(), ofstream::out | ofstream::app);
			time_t rawtime;
		  char time_buffer [80];
			time ( &rawtime );
			strftime (time_buffer,80,"%x %X", localtime ( &rawtime ) );
			if (free==0) outfile << time_buffer << " :: " << toolname_ << " unsuccessfully requested ID (pool is empty!)\n";
			else outfile << time_buffer << " :: " << toolname_ << " requested ID '" << id << "'\n";
			outfile.close();
		}

		// release lock file
		flock.unlock();

		return true;
	}

	bool DocumentIDTagger::tag(DocumentIdentifier& map) const
	{
		String id="";Int free(0);
		try
		{
			if (getID_(id, free, false))
			{	
				if (free>0)
				{
					map.setIdentifier(id);
					return true;
				}
			}
		}
		catch(...) {}

		map.setIdentifier("InvalidID");

		String msg;
		if (free==0) msg = String("Tool ")+toolname_+String(" requested identifier from depleted ID pool '") + getPoolFile() + String("'");
		else msg = String("Tool ")+toolname_+String(" requested identifier from unaccessible ID pool '") + getPoolFile() + String("'. There should be ") + String(free) + String(" identifiers available!");

		throw Exception::DepletedIDPool(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IDTagger", msg);
	}

	bool DocumentIDTagger::countFreeIDs(Int& free) const
	{
		String id="";
		try
		{
			if (getID_(id, free, true))
			{		
				return true;
			}
		}
		catch(...) {}

		return false;
	}

} // namespace OpenMS

