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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>

#include <OpenMS/CONCEPT/StreamHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

using std::ios_base;
using std::ostringstream;
using std::ofstream;

namespace OpenMS
{
  StreamHandler::StreamHandler()
  {

  }

  StreamHandler::StreamHandler(const StreamHandler & source)
  {
    name_to_stream_map_ = source.name_to_stream_map_;
    name_to_counter_map_ = source.name_to_counter_map_;
    name_to_type_map_ = source.name_to_type_map_;
  }

  StreamHandler::~StreamHandler()
  {
    // close all associated streams
    for (map<String, ostream *>::iterator iter = name_to_stream_map_.begin(); iter != name_to_stream_map_.end(); ++iter)
    {
      ostream * stream_pointer = iter->second;
      // file streams need to be closed before
      if (name_to_type_map_[iter->first] == FILE)
      {
        (static_cast<ofstream *>(stream_pointer))->close();
      }
      delete stream_pointer; // call destructor
    }
  }

  StreamHandler & StreamHandler::operator=(const StreamHandler & source)
  {
    name_to_stream_map_ = source.name_to_stream_map_;
    name_to_counter_map_ = source.name_to_counter_map_;
    name_to_type_map_ = source.name_to_type_map_;
    return *this;
  }

  ostream & StreamHandler::getStream(StreamType const type, const String & stream_name)
  {
    if (hasStream(type, stream_name))
    {
      return *name_to_stream_map_[stream_name];
    }
    else
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream_name);
    }
  }

  ostream * StreamHandler::createStream_(const StreamType type, const String & stream_name)
  {
    ostream * stream_pointer;
    switch (type)
    {
    case STRING:
      stream_pointer = new ostringstream();
      break;

    case FILE:
    default:
      stream_pointer = new ofstream(File::absolutePath(stream_name).c_str(), ios_base::app);
      break;
    }

    return stream_pointer;
  }

  Int StreamHandler::registerStream(StreamType const type, const String & stream_name)
  {
    Int state = 1;

    if (name_to_stream_map_.count(stream_name) == 0) // this is an unknown stream .. register
    {
      name_to_stream_map_[stream_name] = createStream_(type, stream_name);
      name_to_type_map_[stream_name] = type;
      name_to_counter_map_[stream_name] = 1;

      // check stream
      if (name_to_stream_map_[stream_name]->fail())
      {
        state = 1; // indicate that something went wrong while creating this stream
      }
    }
    else
    {
      // check type consistency
      if (name_to_type_map_[stream_name] != type)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "This stream was already registered with a different type.");
      }
      ++name_to_counter_map_[stream_name];
    }

    return state;
  }

  bool StreamHandler::hasStream(const StreamType type, const String & stream_name)
  {
    if (name_to_stream_map_.count(stream_name) != 0)
    {
      return name_to_type_map_[stream_name] == type;
    }
    else
    {
      return false;
    }
  }

  void StreamHandler::unregisterStream(StreamType const type, const String & stream_name)
  {
    if (name_to_stream_map_.count(stream_name) != 0) // check if we know this stream
    {
      if (name_to_counter_map_[stream_name] > 1)
      {
        // if there are still references left to this stream
        // just decrease the number of references
        --name_to_counter_map_[stream_name];
      }
      else
      {
        // delete the stream
        if (type == FILE)
        {
          // file streams need to be closed before
          (static_cast<ofstream *>(name_to_stream_map_[stream_name]))->close();
        }

        delete name_to_stream_map_[stream_name];

        // remove entry from the local registry
        name_to_stream_map_.erase(stream_name);
        name_to_counter_map_.erase(stream_name);
        name_to_type_map_.erase(stream_name);
      }
    }
    else
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream_name);
    }
  }

  std::ostream & operator<<(std::ostream & os, StreamHandler const & stream_handler)
  {
    for (map<String, ostream *>::const_iterator iter = stream_handler.name_to_stream_map_.begin(); iter != stream_handler.name_to_stream_map_.end(); ++iter)
    {
      os << "[" << iter->first << "] of type";

      if ((stream_handler.name_to_type_map_.find(iter->first))->second == StreamHandler::FILE)
      {
        os << " FILE";
      }
      else
      {
        os << " STRING";
      }
      os << " #" << (stream_handler.name_to_counter_map_.find(iter->first))->second << " " << iter->second << std::endl;
    }
    return os;
  }

} // end namespace OpenMS
