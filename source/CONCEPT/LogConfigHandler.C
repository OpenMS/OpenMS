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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <iostream>

#include <OpenMS/CONCEPT/LogConfigHandler.h>
#include <OpenMS/CONCEPT/StreamHandler.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace OpenMS::Logger;
using std::cout;
using std::cerr;
using std::endl;

namespace OpenMS
{
  String LogConfigHandler::PARAM_NAME = "log";

  LogConfigHandler::LogConfigHandler()
  {
    // add default configuration
    fatal_streams_.insert("cerr");
    error_streams_.insert("cerr");
    
    warn_streams_.insert("cout");
    info_streams_.insert("cout");
  }

  LogConfigHandler::LogConfigHandler(const LogConfigHandler &other)
  {
    debug_streams_ = other.debug_streams_;
    info_streams_ = other.info_streams_;
    warn_streams_ = other.warn_streams_;
    error_streams_ = other.error_streams_;
    fatal_streams_ = other.fatal_streams_;

    stream_type_map_ = other.stream_type_map_;
  }

  LogConfigHandler::~LogConfigHandler()
  {
  }

  LogConfigHandler& LogConfigHandler::operator= (const LogConfigHandler& source)
  {
    debug_streams_ = source.debug_streams_;
    info_streams_ = source.info_streams_;
    warn_streams_ = source.warn_streams_;
    error_streams_ = source.error_streams_;
    fatal_streams_ = source.fatal_streams_;

    stream_type_map_ = source.stream_type_map_;

    return *this;
  }
  
  Param LogConfigHandler::parse(const StringList & settings)
  {
    Param p;    
    String suffix = " FILE";
    StringList commands;
    for(StringList::ConstIterator iter = settings.begin(); iter != settings.end(); ++iter)
    {
      // split by " " to get all keywords
      StringList l;
      (*iter).split(' ', l, true);
      
      if(l.size() < 2 || l.size() > 3)
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, (*iter), "Error while parsing logger config. Setting can only have 2 or 3 arguments.");
      }

      // we parse a command line here, so we append a FILE to each of the arguments
      // to indicate, that all of these streams are FILE streams
      // for cout/cerr the type parameter is ignored
      String new_command = *iter + suffix;
      commands.push_back(new_command);
    }
    
    p.setValue(LogConfigHandler::PARAM_NAME, commands, "List of all settings that should be applied to the current Logging Configuration");
    
    return p;
  }
  
  void LogConfigHandler::configure(const Param & param)
  {
    StringList configurations = (StringList) param.getValue(LogConfigHandler::PARAM_NAME);
    
    for(StringList::ConstIterator iter = configurations.begin() ; iter != configurations.end() ; ++iter)
    {
      // split by " " to get the commands
      StringList commands;
      iter->split(' ', commands, true );
      
      LogStream & log = getLogStreamByName_(commands[0]);
      
      // convenience variables
      String& command = commands[1];

      // identify action
      if(command == "add")
      {
        // convenience variables
        const String& stream_name = commands[2];

        // add the stream given by the 3rd argument to the defined log
        if(stream_name == "cout")
        {          
          log.insert(cout);
        }
        else if(stream_name == "cerr")
        {
          log.insert(cerr);
        }
        else 
        {
					if (commands.size() <= 3)
					{	// write error to cerr and not a LogStream (because we're just configuring it...)
						std::cerr << "Error during configuring logging: the command '" << (*iter) << "' requires 4 entries but has only " << commands.size() << "\n";
						continue;
					}
					const String& stream_type = commands[3];
						
					// check if a stream with the same name, but different type was already registered
					if(stream_type_map_.count(stream_name) != 0)
					{
						if(stream_type_map_[stream_name] != getStreamTypeByName_(stream_type))
						{
							throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "A stream with the same name but different type was already registered.");
						}
					}
					
					StreamHandler::StreamType type = getStreamTypeByName_(stream_type);
          Int status = STREAM_HANDLER.registerStream(type, stream_name);

          if(!status)
          {
            // operation failed
            throw Exception::FileNotWritable(__FILE__, __LINE__, __PRETTY_FUNCTION__, commands[2]);
          }

          log.insert(STREAM_HANDLER.getStream(type, stream_name));
          log.setPrefix(STREAM_HANDLER.getStream(type, stream_name), "[%S] ");

          stream_type_map_[stream_name] = type;

        }

        // register the stream internally, so that the LogConfigHandler knows
        // which streams were added
        getConfigSetByName_(commands[0]).insert(stream_name);
      }
      else if(command == "remove")
      {
        // convenience variables
        const String& stream_name = commands[2];

        // add the stream given by the 3rd argument to the defined log
        if(stream_name == "cout")
        {
          log.remove(cout);
        }
        else if(stream_name == "cerr")
        {
          log.remove(cerr);
        }
        else 
        {
					if (commands.size() <= 3)
					{	// write error to cerr and not a LogStream (because we're just configuring it...)
						std::cerr << "Error during configuring logging: the command '" << (*iter) << "' requires 4 entries but has only " << commands.size() << "\n";
						continue;
					}
					const String& stream_type = commands[3];
          StreamHandler::StreamType type = getStreamTypeByName_(stream_type);

          // it is a file, get the ostream from the StreamHandler
          if(STREAM_HANDLER.hasStream(type, stream_name))
          {
            log.remove(STREAM_HANDLER.getStream(type, stream_name));
            STREAM_HANDLER.unregisterStream(type, stream_name);
          }
        }
        // unregister the stream internally, so that the LogConfigHandler knows
        // which streams were added
        getConfigSetByName_(commands[0]).erase(stream_name);

        // remove the type from the stream_type_map if there is no
        // stream referencing it anymore
        if(! STREAM_HANDLER.hasStream(stream_type_map_[stream_name], stream_name))
        {
          stream_type_map_.erase(stream_name);
        }
      }
      else if(command == "clear")
      {
        // remove all streams from the given log
        for(std::set<String>::iterator it = getConfigSetByName_(commands[0]).begin() ; it != getConfigSetByName_(commands[0]).end() ; ++it)
        {
          if(*it == "cout")
          {
            log.remove(cout);
          }
          else if(*it == "cerr")
          {
            log.remove(cerr);
          }
          else // handle the file streams
          {
            log.remove(STREAM_HANDLER.getStream(stream_type_map_[*it], *it));
            STREAM_HANDLER.unregisterStream(stream_type_map_[*it], *it);

            // remove the type from the stream_type_map if there is no
            // stream referencing it anymore
            if(! STREAM_HANDLER.hasStream(stream_type_map_[*it], *it))
            {
              stream_type_map_.erase(*it);
            }
          }
        }

        // clean the set
        getConfigSetByName_(commands[0]).clear();
      }
    }
  }
    
  LogStream & LogConfigHandler::getLogStreamByName_(const String & stream_name)
  {
    LogStream * log = &Log_debug; // default
    
    if(stream_name == "DEBUG")
    {
      log = &Log_debug;
    }
    else if(stream_name == "INFO")
    {
      log = &Log_info;
    }
    else if(stream_name == "WARNING")
    {
      log = &Log_warn;
    }
    else if(stream_name == "ERROR")
    {
      log = &Log_error;
    }
    else if(stream_name == "FATAL_ERROR")
    {
      log = &Log_fatal;
    }
    else
    {
      Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream_name);
    }
    return *log;
  }
  
  std::set<String> & LogConfigHandler::getConfigSetByName_(const String & stream_type)
  {
    std::set<String> * s = &debug_streams_;
    if(stream_type == "DEBUG")
    {
      s = &debug_streams_;
    }
    else if(stream_type == "INFO")
    {
      s = &info_streams_;
    }
    else if(stream_type == "WARNING")
    {
      s = &warn_streams_;
    }
    else if(stream_type == "ERROR")
    {
      s = &error_streams_;
    }
    else if(stream_type == "FATAL_ERROR")
    {
      s = &fatal_streams_;
    }
    else
    {
      Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream_type);
    }  
    
    return *s;
  }

  ostream& LogConfigHandler::getStream(const String &name)
  {
    if(stream_type_map_.count(name) != 0)
    {
      return STREAM_HANDLER.getStream(stream_type_map_[name], name);
    }
    else
    {
      // there is no stream with this name
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "There is no stream with the given name.");
    }
  }

  StreamHandler::StreamType LogConfigHandler::getStreamTypeByName_(const String &stream_type)
  {
    StreamHandler::StreamType type;
    if(stream_type == "FILE")
    {
      type = StreamHandler::FILE;
    }
    else if(stream_type == "STRING")
    {
      type = StreamHandler::STRING;
    }
    else
    {
      // unsupported log type
      throw Exception::IllegalArgument(__FILE__, __LINE__ , __PRETTY_FUNCTION__ , "The log type " + stream_type + " is not supported");
    }

    return type;
  }
  
  void printStreamConfig_(std::ostream& os, const String & name, const std::set<String> & stream_names, const std::map<String,StreamHandler::StreamType> & stream_type_map)
  {
    os << name << endl;
    for(std::set<String>::const_iterator it = stream_names.begin() ; it != stream_names.end() ; ++it)
    {
      os << "->" << "\t" << *it;
      // append stream type
      os << " (";

      switch((stream_type_map.find(*it))->second)
      {
      case StreamHandler::STRING:
        os << "STRINGSTREAM";
        break;
      case StreamHandler::FILE:
      default:
        os << "FILE";
        break;
      }
      os <<  ")";
      os << std::endl;
    }
  }
  
  std::ostream& operator<<(std::ostream& os, LogConfigHandler const & lch)
  {

    printStreamConfig_(os, "LOG_DEBUG", lch.debug_streams_, lch.stream_type_map_);
    printStreamConfig_(os, "LOG_INFO", lch.info_streams_, lch.stream_type_map_);
    printStreamConfig_(os, "LOG_WARNING", lch.warn_streams_, lch.stream_type_map_);
    printStreamConfig_(os, "LOG_ERROR", lch.error_streams_, lch.stream_type_map_);
    printStreamConfig_(os, "LOG_FATAL_ERROR", lch.fatal_streams_, lch.stream_type_map_);
    
    return os;
  }

  LogConfigHandler* LogConfigHandler::instance_ = NULL;

  LogConfigHandler& LogConfigHandler::getInstance()
  {
    if(LogConfigHandler::instance_ == 0)
    {
      LogConfigHandler::instance_ = new LogConfigHandler();
    }
    return *LogConfigHandler::instance_;
  }
} // end namespace OpenMS
