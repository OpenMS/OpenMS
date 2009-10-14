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
// $Maintainer: Chris Bielow, Stephan Aiche, Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------


/**

	Generously provided by the BALL people, taken from version 1.2

	Originally implemented by OK who refused to take any responsibility
	for the code ;)
*/
#include <limits>
#include <string>
#include <cstring>
#include <cstdio>
#include <OpenMS/CONCEPT/LogStream.h>

#define BUFFER_LENGTH 32768

using namespace std;

namespace OpenMS 
{
	namespace Logger
	{

	const LogLevel LogStreamBuf::MIN_LEVEL = DEVELOPMENT;
	const LogLevel LogStreamBuf::MAX_LEVEL = FATAL_ERROR;
	const time_t LogStreamBuf::MAX_TIME = numeric_limits<time_t>::max();

	LogStreamBuf::LogStreamBuf() 
		: std::streambuf(),
			pbuf_(0),
			level_(DEVELOPMENT),
			tmp_level_(DEVELOPMENT),
			stream_list_(),
			incomplete_line_(),
      log_cache_counter_(0),
      log_cache_(),
      log_time_cache_()
	{
		pbuf_ = new char [BUFFER_LENGTH];
		std::streambuf::setp(pbuf_, pbuf_ + BUFFER_LENGTH - 1);
	}
		
	LogStreamBuf::~LogStreamBuf() 
	{
		sync();
    clearCache_();
		delete [] pbuf_;
	}


  int LogStreamBuf::overflow(int c)
  {
    sync();
    return ::std::streambuf::overflow(c);
  }
		

  LogStreamBuf* LogStream::rdbuf() 
  {
    return (LogStreamBuf*)std::ios::rdbuf();
  }


  LogStreamBuf* LogStream::operator -> () 
  {
    return rdbuf();
  }


  void LogStream::setLevel(LogLevel level) 
  {
    if (rdbuf() == 0)
    {
      return;
    }

    // set the new level
    rdbuf()->level_ = level;
  
    // set tmp_level_, too - to otherwise the
    // new level would take effect in the line after 
    // the next!
    rdbuf()->tmp_level_ = level;
  }


  LogLevel LogStream::getLevel() 
  {
    if (rdbuf() != 0)
    {
      return rdbuf()->level_;
    }
    else 
    {
      return DEVELOPMENT;
    }
  }


  LogStream& LogStream::level(LogLevel level) 
  {
    // set the temporary level 
    // will be reset by sync(), i.e. at the end of the next line
    if (rdbuf() != 0)
    {
      rdbuf()->tmp_level_ = level;
    }
    
    return *this;
  }


  // caching methods
  Size LogStreamBuf::getNextLogCounter_()
  {
    return ++log_cache_counter_;
  }


  bool LogStreamBuf::isInCache_(std::string const & line)
  {
    if(log_cache_.count(line) == 0) return false;
    else {
      // increment counter
      log_cache_[line].counter++;
      
      // remove old entry 
      log_time_cache_.erase(log_cache_[line].timestamp);

      // update timestamp
      Size counter_value = getNextLogCounter_();
      log_cache_[line].timestamp = counter_value;
      log_time_cache_[counter_value] = line;
      return true;
    }
  }
 
  std::string LogStreamBuf::addToCache_(std::string const & line)
  {
    std::string extra_message = "";

    if(log_cache_.size() > 1) // check if we need to remove one of the entries
    {
      // get smallest key
      map<Size, string>::iterator it = log_time_cache_.begin();

      // check if message occured more then once
      if(log_cache_[it->second].counter != 0) {
        std::stringstream stream;
        stream << "<" << it->second << "> occured " << ++log_cache_[it->second].counter << " times";
        extra_message = stream.str();
      }

      log_cache_.erase(it->second);
      log_time_cache_.erase(it);
    }

    Size counter_value =getNextLogCounter_();
    log_cache_[line].counter = 0;
    log_cache_[line].timestamp = counter_value;

    log_time_cache_[counter_value] = line;

    return extra_message;
  }


  void LogStreamBuf::clearCache_() 
  {
    // if there are any streams in our list, we
    // copy the line into that streams, too and flush them
    map<std::string, LogCacheStruct >::iterator it = log_cache_.begin();

    for(; it != log_cache_.end() ; ++it) 
    {
      if((it->second).counter != 0) 
      {
        std::stringstream stream;
        stream << "<" << it->first << "> occured " << ++(it->second).counter << " times";
        std::list<StreamStruct>::iterator list_it = stream_list_.begin();
        for (; list_it != stream_list_.end(); ++list_it)
        {
          // if the stream is open for that level, write to it...
          if ((list_it->min_level <= tmp_level_) && (list_it->max_level >= tmp_level_))
          {
            *(list_it->stream) << expandPrefix_(list_it->prefix, tmp_level_, time(0)).c_str()
                               << stream.str().c_str() << std::endl;

            if (list_it->target != 0)
            {
              list_it->target->logNotify();
            }
          }
        }
      }
    }
  }

	int LogStreamBuf::sync() 
	{
		static char buf[BUFFER_LENGTH];

		// sync our streambuffer...
		if (pptr() != pbase()) 
		{
				
			char*	line_start = pbase();
			char*	line_end = pbase();

			while (line_end < pptr())
			{
				// search for the first end of line
				for (; line_end < pptr() && *line_end != '\n'; line_end++) {};

				if (line_end >= pptr()) 
				{
					// Copy the incomplete line to the incomplete_line_ buffer
					size_t length = line_end - line_start + 1;
					length = std::min(length, (size_t)(BUFFER_LENGTH - 1));
					strncpy(&(buf[0]), line_start, length);

					// if length was too large, we copied one byte less than BUFFER_LENGTH to have
					// room for the final \0
					buf[length] = '\0';

					incomplete_line_ += &(buf[0]);

					// mark everything as read
					line_end = pptr() + 1;
				} 
				else 
				{
					// note: pptr() - pbase() should be bounded by BUFFER_LENGTH, so this should always work
					memcpy(&(buf[0]), line_start, line_end - line_start + 1);
					buf[line_end - line_start] = '\0';
						
					// assemble the string to be written
					// (consider leftovers of the last buffer from incomplete_line_)
					std::string outstring = incomplete_line_;
					incomplete_line_ = "";
					outstring += &(buf[0]);

          // chech if we already have that in line in our cache
          if(!isInCache_(outstring)) 
          {

            // add line to the log cache
            std::string extra_message = addToCache_(outstring);

            // if there are any streams in our list, we
            // copy the line into that streams, too and flush them
            std::list<StreamStruct>::iterator list_it = stream_list_.begin();
            for (; list_it != stream_list_.end(); ++list_it)
            {
              // if the stream is open for that level, write to it...
              if ((list_it->min_level <= tmp_level_) && (list_it->max_level >= tmp_level_))
              {
                // first print the message from cache 
                if(extra_message != "") 
                {
                  *(list_it->stream) << expandPrefix_(list_it->prefix, tmp_level_, time(0)).c_str()
                                     << extra_message.c_str() << std::endl;
                }

                *(list_it->stream) << expandPrefix_(list_it->prefix, tmp_level_, time(0)).c_str()
                                   << outstring.c_str() << std::endl;

                if (list_it->target != 0)
                {
                  list_it->target->logNotify();
                }
              }
            }
			
            // update the line pointers (increment both)
            line_start = ++line_end;
					
            // remove cr/lf from the end of the line				
            while (outstring.size() && (outstring[outstring.size() - 1] == 10 || outstring[outstring.size() - 1] == 13))
            {
              std::string::iterator p = outstring.end();
              p--;
              outstring.erase(p);
            }
		
            // reset tmp_level_ to the previous level
            // (needed for LogStream::level() only)
            tmp_level_ = level_;
          } else { line_start = ++line_end; }
				}
			}

			// remove all processed lines from the buffer
			pbump((int)(pbase() - pptr()));
		}
	
		return 0;
	}

	string LogStreamBuf::expandPrefix_
		(const std::string& prefix, LogLevel level, time_t time) const
	{
		string::size_type	index = 0;
		Size copied_index = 0;
		string result("");

		while ((index = prefix.find("%", index)) != string::npos) 
		{
			// append any constant parts of the string to the result
			if (copied_index < index) 
			{
				result.append(prefix.substr(copied_index, index - copied_index));
				copied_index = (SignedSize)index;
			}
			
			if (index < prefix.size()) 
			{
				char	buffer[64];
				char*	buf = &(buffer[0]);

				switch (prefix[index + 1]) 
				{
					case '%': // append a '%' (escape sequence)
						result.append("%");
						break;

					case 'l': // append the loglevel
						sprintf(buf, "%d", level);
						result.append(buf);
						break;

					case 'y':	// append the message type (error/warning/information)
						result.append(LogLevelToStringUpper(level));
						break;

					case 'T':	// time: HH:MM:SS
						strftime(buf, BUFFER_LENGTH - 1, "%H:%M:%S", localtime(&time));
						result.append(buf);
						break;

					case 't': // time: HH:MM	
						strftime(buf, BUFFER_LENGTH - 1, "%H:%M", localtime(&time));
						result.append(buf);
						break;

					case 'D':	// date: DD.MM.YYYY
						strftime(buf, BUFFER_LENGTH - 1, "%d.%m.%Y", localtime(&time));
						result.append(buf);
						break;

					case 'd':	// date: DD.MM.
						strftime(buf, BUFFER_LENGTH - 1, "%d.%m.", localtime(&time));
						result.append(buf);
						break;

					case 'S':	// time+date: DD.MM.YYYY, HH:MM:SS
						strftime(buf, BUFFER_LENGTH - 1, "%d.%m.%Y, %H:%M:%S", localtime(&time));
						result.append(buf);
						break;

					case 's':	// time+date: DD.MM., HH:MM
						strftime(buf, BUFFER_LENGTH - 1, "%d.%m., %H:%M", localtime(&time));
						result.append(buf);
						break;

					default:
						break;
				}
				index += 2;
				copied_index += 2;
			}
		}

		if (copied_index < prefix.size()) 
		{
			result.append(prefix.substr(copied_index, prefix.size() - copied_index));
		}

		return result;
	}

	LogStreamNotifier::LogStreamNotifier()
		: registered_at_(0)
	{
	}

	LogStreamNotifier::~LogStreamNotifier()
	{
		unregister();
	}

	void LogStreamNotifier::logNotify()
	{
	}

	void LogStreamNotifier::unregister()
	{
		if (registered_at_ == 0) return;

		registered_at_->remove(stream_);
		registered_at_ = 0;
	}

	void LogStreamNotifier::registerAt(LogStream& log, LogLevel min_level, LogLevel max_level)
	{
		unregister();

		registered_at_ = &log;
		log.insertNotification(stream_, *this, min_level, max_level);
	}

	// keep the given buffer	
	LogStream::LogStream(LogStreamBuf* buf, bool delete_buf, bool associate_stdio)
		: std::ios(buf),
			std::ostream(buf),
			delete_buffer_(delete_buf),
			disable_output_(false)
	{
		if (associate_stdio) 
		{
			// associate cout to informations and warnings,
			// cerr to errors by default
			insert(std::cout, WARNING, FATAL_ERROR);
			insert(std::cerr, ERROR);
		}
	}

	LogStream::~LogStream()
	{
		if (delete_buffer_)
		{
			// remove the streambuffer
			delete rdbuf();
		}
	}

	void LogStream::insert(std::ostream& stream, LogLevel min_level, LogLevel max_level) 
	{
		if (!bound_() || hasStream_(stream))
		{
			return;
		}
			
		// we didn`t find it - create a new entry in the list
		LogStreamBuf::StreamStruct s_struct;
		s_struct.min_level = min_level;
		s_struct.max_level = max_level;
		s_struct.stream = &stream;
		rdbuf()->stream_list_.push_back(s_struct);
	}

	void LogStream::remove(std::ostream& stream) 
	{
		if (!bound_()) return;

		StreamIterator it = findStream_(stream);
		if (it != rdbuf()->stream_list_.end())
		{
			rdbuf()->stream_list_.erase(it);
		}
	}

	void LogStream::insertNotification(std::ostream& s, LogStreamNotifier& target,
																		LogLevel min_level, LogLevel max_level)
	{
		if (!bound_()) return;

		insert(s, min_level, max_level);

		StreamIterator it = findStream_(s);
		(*it).target = &target;
	}

	LogStream::StreamIterator LogStream::findStream_(const std::ostream& s)
	{
		StreamIterator list_it = rdbuf()->stream_list_.begin();
		for (; list_it != rdbuf()->stream_list_.end(); ++list_it)
		{
			if (list_it->stream == &s) 
			{
				return list_it;
			}
		}

		return list_it;
	}

	bool LogStream::hasStream_(std::ostream& stream)
	{
		if (!bound_()) return false;

		return findStream_(stream) != rdbuf()->stream_list_.end();
	}

	void LogStream::setMinLevel(const std::ostream& stream, LogLevel level) 
	{
		if (!bound_()) return;
			
		StreamIterator it = findStream_(stream);
		if (it != rdbuf()->stream_list_.end())
		{
			(*it).min_level = level;
		}
	}

	void LogStream::setMaxLevel(const std::ostream& stream, LogLevel level) 
	{
		if (!bound_()) return;
			
		StreamIterator it = findStream_(stream);
		if (it != rdbuf()->stream_list_.end())
		{
			(*it).max_level = level;
		}
	}
	
	void LogStream::setPrefix(const std::ostream& s, const string& prefix) 
	{
		if (!bound_()) return;

		StreamIterator it = findStream_(s);
		if (it != rdbuf()->stream_list_.end())
		{
			(*it).prefix = prefix;
		}		
	}

	void LogStream::setPrefix(const string& prefix)
	{
		if (!bound_()) return;
		for (StreamIterator it = rdbuf()->stream_list_.begin(); it != rdbuf()->stream_list_.end(); ++it)
		{
			(*it).prefix = prefix;
		}
	}
	
	void LogStream::disableOutput()
		
	{
		disable_output_ = true;
	}

	void LogStream::enableOutput()
		
	{
		disable_output_ = false;
		std::ostream::flush();
	}

	bool LogStream::outputEnabled() const
		
	{
		return disable_output_;
	}

	void LogStream::flush()
		
	{
		if (disable_output_) return;

		std::ostream::flush();
	}

	bool LogStream::bound_() const
	{
		LogStream*	non_const_this = const_cast<LogStream*>(this);

		return (non_const_this->rdbuf() != 0);
	}
	
	} // namespace Logger

	// global default logstream
	OPENMS_DLLAPI	Logger::LogStream Log(new Logger::LogStreamBuf, true, true);

} // namespace OpenMS
