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

#include <cstdio>
#include <limits>
#include <cstring>
#include <OpenMS/CONCEPT/LogStream.h>

#define BUFFER_LENGTH 32768

using namespace std;

namespace OpenMS 
{
	const Int LogStreamBuf::MIN_LEVEL = numeric_limits<Int>::min();
	const Int LogStreamBuf::MAX_LEVEL = numeric_limits<Int>::max();
	const Time LogStreamBuf::MAX_TIME = numeric_limits<Time>::max();

	LogStreamBuf::LogStreamBuf() 
		: std::streambuf(),
			pbuf_(0),
			loglines_(),
			level_(0),
			tmp_level_(0),
			stream_list_(),
			incomplete_line_()
	{
		pbuf_ = new char [BUFFER_LENGTH];
		std::streambuf::setp(pbuf_, pbuf_ + BUFFER_LENGTH - 1);
	}
		
	LogStreamBuf::~LogStreamBuf() 
	{
		sync();

		delete [] pbuf_;
	}

	void LogStreamBuf::dump(std::ostream& stream) 
	{
		char buf[BUFFER_LENGTH];
		Size line;
		for (line = (Size)loglines_.size(); line > 0; --line) 
		{
			strftime(&(buf[0]), BUFFER_LENGTH - 1, "%d.%m.%Y %H:%M:%S ", localtime(&(loglines_[line - 1].time)));
			stream << buf << "[" << loglines_[line - 1].level
						 << "]:" << loglines_[line - 1].text.c_str() << std::endl;
		}
	}
 
	Int LogStreamBuf::sync() 
	{
		static char buf[BUFFER_LENGTH];

		// sync our streambuffer...
		if (pptr() != pbase()) 
		{
				
			char*	line_start = pbase();
			char*	line_end = pbase();
			
			while (line_end <= pptr())
			{
				// search for the first end of line
				for (; line_end < pptr() && *line_end != '\n'; line_end++) ;

				if (line_end >= pptr()) 
				{
					// Copy the incomplete line to the incomplete_line_ buffer
					size_t length = line_end - line_start + 1;
					length = std::max(length, (size_t)(BUFFER_LENGTH - 1));
					strncpy(&(buf[0]), line_start, length);
					buf[line_end - line_start] = '\0';
					incomplete_line_ += &(buf[0]);

					// mark everything as read
					line_end = pptr() + 1;
				} 
				else 
				{
					memcpy(&(buf[0]), line_start, line_end - line_start + 1);
					buf[line_end - line_start] = '\0';
						
					// assemble the String to be written
					// (consider leftovers of the last buffer from incomplete_line_)
					String outstring = incomplete_line_;
					incomplete_line_ = "";
					outstring += &(buf[0]);

					// if there are any streams in our list, we
					// copy the line into that streams, too and flush them
					std::list<StreamStruct>::iterator list_it = stream_list_.begin();
					for (; list_it != stream_list_.end(); ++list_it)
					{
						// if the stream is open for that level, write to it...
						if ((list_it->min_level <= tmp_level_) && (list_it->max_level >= tmp_level_))
						{
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
						String::iterator p = outstring.end();
						p--;
						outstring.erase(p);
					}
		
					// store the line 
					Logline	logline;

					logline.text = outstring;
					logline.level = tmp_level_;
					logline.time = time(0);
			
					// store the new line
					loglines_.push_back(logline);

					// reset tmp_level_ to the previous level
					// (needed for LogStream::level() only)
					tmp_level_ = level_;
				}
			}

			// remove all processed lines from the buffer
			pbump((Int)(pbase() - pptr()));
		}
	
		return 0;
	}

	String LogStreamBuf::expandPrefix_
		(const String& prefix, Int level, Time time) const
	{
		String::size_type	index = 0;
		Size copied_index = 0;
		String result("");

		while ((index = prefix.find("%", index)) != String::npos) 
		{
			// append any constant parts of the String to the result
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
						if (level >= LogStream::ERROR_LEVEL) 
						{
							result.append("ERROR");
						}
						else 
						{
							if (level >= LogStream::WARNING_LEVEL) 
							{
								result.append("WARNING");
							}
							else 
							{
								if (level >= LogStream::INFORMATION_LEVEL) 
								{
									result.append("INFORMATION");
								}
								else 
								{
									result.append("LOG");
								}
							}
						}
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

	void LogStreamNotifier::registerAt(LogStream& log, Int min_level, Int max_level)
	{
		unregister();

		registered_at_ = &log;
		log.insertNotification(stream_, *this, min_level, max_level);
	}

	// keep the given buffer	
	LogStream::LogStream(LogStreamBuf* buf, bool delete_buf, bool associate_stdio)
		: basic_ios<char>(buf),
			ostream(buf),
			delete_buffer_(delete_buf),
			disable_output_(false)
	{
		if (associate_stdio) 
		{
			// associate cout to informations and warnings,
			// cerr to errors by default
			insert(std::cout, INFORMATION_LEVEL, ERROR_LEVEL - 1);
			insert(std::cerr, ERROR_LEVEL);
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
	
	void LogStream::clear()
	{
		rdbuf()->loglines_.clear();
	}

	void LogStream::insert(std::ostream& stream, Int min_level, Int max_level) 
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
																		Int min_level, Int max_level)
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

	void LogStream::setMinLevel(const std::ostream& stream, Int level) 
	{
		if (!bound_()) return;
			
		StreamIterator it = findStream_(stream);
		if (it != rdbuf()->stream_list_.end())
		{
			(*it).min_level = level;
		}
	}

	void LogStream::setMaxLevel(const std::ostream& stream, Int level) 
	{
		if (!bound_()) return;
			
		StreamIterator it = findStream_(stream);
		if (it != rdbuf()->stream_list_.end())
		{
			(*it).max_level = level;
		}
	}
	
	void LogStream::setPrefix(const std::ostream& s, const String& prefix) 
	{
		if (!bound_()) return;

		StreamIterator it = findStream_(s);
		if (it != rdbuf()->stream_list_.end())
		{
			(*it).prefix = prefix;
		}		
	}
	
	Size LogStream::getNumberOfLines(Int min_level, Int max_level) const  
	{
		if (!bound_()) return 0;

		// iterate over all loglines and count the lines of interest
		LogStream*	non_const_this = const_cast<LogStream*>(this);
		vector<LogStreamBuf::Logline>::iterator	it = non_const_this->rdbuf()->loglines_.begin();
		Size	count = 0;

		for (; it != non_const_this->rdbuf()->loglines_.end(); ++it) 
		{
			if ((*it).level >= min_level && (*it).level <= max_level)
			{
				count++;
			}
		}

		return count;
	}

	String LogStream::getLineText(const SignedSize& index) const
	{
		if ((SignedSize)getNumberOfLines() < index)
		{
			return "";
		}

		if (!bound_()) return "";

		return const_cast<LogStream*>(this)->rdbuf()->loglines_[index].text;	
	}

	Int LogStream::getLineLevel(const SignedSize& index) const
	{
		if ((SignedSize)getNumberOfLines() < index)
		{
			return -1;
		}

		if (!bound_()) return -1;

		return const_cast<LogStream*>(this)->rdbuf()->loglines_[index].level;	
	}


	time_t LogStream::getLineTime(const SignedSize& index) const
	{
		if ((SignedSize)getNumberOfLines() < index)
		{
			return 0;
		}

		if (!bound_()) return 0;

		return const_cast<LogStream*>(this)->rdbuf()->loglines_[index].time;	
	}

	std::list<Int>	LogStream::filterLines
		(Int min_level, Int max_level,
		 Time earliest, Time latest, const String& s) const
	{
		std::list<Int>	list_indices;
		Size pos = 0;
		LogStreamBuf* log = const_cast<LogStream*>(this)->rdbuf();

		while (pos < log->loglines_.size() && 
					 log->loglines_[pos].time < earliest)
		{
			pos++;
		}
		while (pos < log->loglines_.size() && 
					 log->loglines_[pos].time <= latest)
		{
			if (log->loglines_[pos].level >= min_level &&
					log->loglines_[pos].level <= max_level)
			{
				if (s.length() > 0)
				{
					if (log->loglines_[pos].text.find(s, 0) != String::npos)
					{
						list_indices.push_back((Int)pos);
					}
				}
				else
				{
					list_indices.push_back((Int)pos);
				}
			}
			pos++;
		}
		return list_indices;
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

	// global default logstream
	//OPENMS_DLLAPI	LogStream	Log(new LogStreamBuf, true, true);

  int LogStreamBuf::overflow(Int c)
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

  void LogStream::setLevel(Int level) 
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

  Int LogStream::getLevel() 
  {
    if (rdbuf() != 0)
    {
      return rdbuf()->level_;
    }
    else 
    {
      return 0;
    }
  }

  LogStream& LogStream::level(Int level) 
  {
    // set the temporary level 
    // will be reset by sync(), i.e. at the end of the next line
    if (rdbuf() != 0)
    {
      rdbuf()->tmp_level_ = level;
    }

    return *this;
  }

  LogStream& LogStream::error(Int level)
  {
    // set the temporary level to ERROR
    // will be reset by sync(), i.e. at the end of the next line
    if (rdbuf() != 0)
    {
      rdbuf()->tmp_level_ = ERROR_LEVEL + level;
    }

    return *this;
  }

  LogStream& LogStream::warn(Int level)
  {
    // set the temporary level to WARNING
    // will be reset by sync(), i.e. at the end of the next line
    if (rdbuf() != 0)
    {
      rdbuf()->tmp_level_ = WARNING_LEVEL + level;
    }

    return *this;
  }

  LogStream& LogStream::info(Int level)
  {
    // set the temporary level to INFORMATION
    // will be reset by sync(), i.e. at the end of the next line
    if (rdbuf() != 0)
    {
      rdbuf()->tmp_level_ = INFORMATION_LEVEL + level;
    }

    return *this;
  }

}
