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
  
  StreamHandler::StreamHandler(const StreamHandler& source)
  {
    name_to_stream_map_ = source.name_to_stream_map_;
    name_to_counter_map_ = source.name_to_counter_map_;
    name_to_type_map_ = source.name_to_type_map_;
  }
  
  StreamHandler::~StreamHandler()
  {
    // close all associated streams
    for(map<String, ostream* >::iterator iter = name_to_stream_map_.begin(); iter != name_to_stream_map_.end() ; ++iter)
    {
      ostream * stream_pointer = iter->second;
      // file streams need to be closed before
      if(name_to_type_map_[iter->first] == FILE)
      {
        (static_cast<ofstream*>(stream_pointer))->close();
      }      
      delete stream_pointer; // call destructor
    }
  }
  
  StreamHandler& StreamHandler::operator =(const StreamHandler& source)
  {
    name_to_stream_map_ = source.name_to_stream_map_;
    name_to_counter_map_ = source.name_to_counter_map_;
    name_to_type_map_ = source.name_to_type_map_;
    return *this;
  }
  
  ostream& StreamHandler::getStream(StreamType const type, const String& stream_name)
  {
    if(hasStream(type, stream_name))
    {
      return *name_to_stream_map_[stream_name];
    }
    else
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream_name);
    }
  }

  ostream* StreamHandler::createStream_(const StreamType type, const String & stream_name)
  {
    ostream* stream_pointer;
    switch(type)
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

    if(name_to_stream_map_.count(stream_name) == 0) // this is an unknown stream .. register
    {
      name_to_stream_map_[stream_name] = createStream_(type, stream_name);
      name_to_type_map_[stream_name] = type;
      name_to_counter_map_[stream_name] = 1;

      // check stream
      if(name_to_stream_map_[stream_name]->fail())
      {
        state = 1; // indicate that something went wrong while creating this stream
      }
    }
    else
    {
      // check type consistency
      if(name_to_type_map_[stream_name] != type)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "This stream was already registered with a different type.");
      }
      ++name_to_counter_map_[stream_name];
    }

    return state;
  }
  
  bool StreamHandler::hasStream(const StreamType type, const String& stream_name)
  {
    if(name_to_stream_map_.count(stream_name) != 0)
    {
      return (name_to_type_map_[stream_name] == type);
    }
    else
    {
      return false;
    }
  }
  
  void StreamHandler::unregisterStream(StreamType const type, const String& stream_name)
  {
    if(name_to_stream_map_.count(stream_name) != 0) // check if we know this stream
    {
      if(name_to_counter_map_[stream_name] > 1)
      {
        // if there are still references left to this stream
        // just decrease the number of references
        --name_to_counter_map_[stream_name];
      }
      else
      {
        // delete the stream
        if(type == FILE)
        {
          // file streams need to be closed before
          (static_cast<ofstream*>(name_to_stream_map_[stream_name]))->close();
        }

        delete name_to_stream_map_[stream_name];

        // remove entry from the local registry
        name_to_stream_map_.erase(stream_name);
        name_to_counter_map_.erase(stream_name);
        name_to_type_map_.erase(stream_name );
      }
    }
    else
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream_name);
    }
  }

  std::ostream& operator<<(std::ostream& os, StreamHandler const & stream_handler)
  {
    for(map<String, ostream* >::const_iterator iter = stream_handler.name_to_stream_map_.begin(); iter != stream_handler.name_to_stream_map_.end() ; ++iter)
    {
      os << "[" << iter->first << "] of type";

      if((stream_handler.name_to_type_map_.find(iter->first))->second == StreamHandler::FILE)
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
