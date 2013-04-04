// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/DocumentIDTagger.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDir>

#ifdef _MSC_VER // disable some boost warnings that distract from ours
#   pragma warning( push ) // save warning state
#   pragma warning( disable : 4018 )
#endif
#include <boost/interprocess/sync/file_lock.hpp>
#ifdef _MSC_VER
#   pragma warning( pop )  // restore old warning state
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
  DocumentIDTagger::DocumentIDTagger(String toolname) :
    toolname_(toolname),
    pool_file_()
  {
    pool_file_ = File::getOpenMSDataPath() + ("/IDPool/IDPool.txt");
  }

  DocumentIDTagger::DocumentIDTagger(const DocumentIDTagger & source) :
    toolname_(source.toolname_),
    pool_file_(source.pool_file_)
  {
  }

  DocumentIDTagger::~DocumentIDTagger()
  {
  }

  DocumentIDTagger & DocumentIDTagger::operator=(const DocumentIDTagger & source)
  {
    if (source == *this)
      return *this;

    toolname_ = source.toolname_;
    pool_file_ = source.pool_file_;
    return *this;
  }

  bool DocumentIDTagger::operator==(const DocumentIDTagger & rhs) const
  {
    return (toolname_ == rhs.toolname_)
           && (pool_file_ == rhs.pool_file_);
  }

  bool DocumentIDTagger::operator!=(const DocumentIDTagger & rhs) const
  {
    return !(operator==(rhs));
  }

  String DocumentIDTagger::getPoolFile() const
  {
    return pool_file_;
  }

  void DocumentIDTagger::setPoolFile(const String & file)
  {
    pool_file_ = file;
  }

  bool DocumentIDTagger::getID_(String & id, Int & free, bool idcount_only) const
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
    if (!in.is_open())
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
    try
    {
      flock.lock();
    }
    catch (exception /*e*/)
    {
      return false;
    }
    // we have the lock!

    // now open temp output file
    ofstream out;
    if (!idcount_only)
    {
      out.open(IDPool_file_tmp.c_str(), ios::out | ios::trunc);
      if (!out.is_open())
      {
        std::cerr << "IDTagger::getID_() " << IDPool_file_tmp << " file failed to open for writing.\n";
        flock.unlock();
        in.close();
        return false;
      }
    }

    // copy the ID pool, but excise the first ID
    string line;
    while (!in.eof())
    {
      getline(in, line);
      if (line.length() == 0)
        continue;
      ++free;
      if (free == 1)
        id = line;                 // pull out first ID
      if (!idcount_only)
      {
        if (free != 1)       // delete first line
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
      rename(IDPool_file_tmp.c_str(), IDPool_file.c_str());

      // write ID to log file
      String logging_file = String(IDPool_file + String(".log"));
      ofstream outfile;
      outfile.open(logging_file.c_str(), ofstream::out | ofstream::app);
      time_t rawtime;
      char time_buffer[80];
      time(&rawtime);
      strftime(time_buffer, 80, "%x %X", localtime(&rawtime));
      if (free == 0)
        outfile << time_buffer << " :: " << toolname_ << " unsuccessfully requested ID (pool is empty!)\n";
      else
        outfile << time_buffer << " :: " << toolname_ << " requested ID '" << id << "'\n";
      outfile.close();
    }

    // release lock file
    flock.unlock();

    return true;
  }

  bool DocumentIDTagger::tag(DocumentIdentifier & map) const
  {
    String id = ""; Int free(0);
    try
    {
      if (getID_(id, free, false))
      {
        if (free > 0)
        {
          map.setIdentifier(id);
          return true;
        }
      }
    }
    catch (...)
    {
    }

    map.setIdentifier("InvalidID");

    String msg;
    if (free == 0)
      msg = String("Tool ") + toolname_ + String(" requested identifier from depleted ID pool '") + getPoolFile() + String("'");
    else
      msg = String("Tool ") + toolname_ + String(" requested identifier from unaccessible ID pool '") + getPoolFile() + String("'. There should be ") + String(free) + String(" identifiers available!");

    throw Exception::DepletedIDPool(__FILE__, __LINE__, __PRETTY_FUNCTION__, "IDTagger", msg);
  }

  bool DocumentIDTagger::countFreeIDs(Int & free) const
  {
    String id = "";
    try
    {
      if (getID_(id, free, true))
      {
        return true;
      }
    }
    catch (...)
    {
    }

    return false;
  }

} // namespace OpenMS
