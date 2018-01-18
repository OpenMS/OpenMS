// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Chris Bielow, Stephan Aiche $
// $Authors: Chris Bielow, Stephan Aiche, Andreas Bertsch $
// --------------------------------------------------------------------------


/**

  Generously provided by the BALL people, taken from version 1.2
  with slight modifications

  Originally implemented by OK who refused to take any responsibility
  for the code ;)
*/
#include <limits>
#include <string>
#include <cstring>
#include <cstdio>
#include <algorithm>    // std::min
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/StreamHandler.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define BUFFER_LENGTH 32768

using namespace std;

namespace OpenMS
{
  namespace Logger
  {

    const time_t LogStreamBuf::MAX_TIME = numeric_limits<time_t>::max();
    const std::string LogStreamBuf::UNKNOWN_LOG_LEVEL = "UNKNOWN_LOG_LEVEL";

    LogStreamBuf::LogStreamBuf(std::string log_level) :
      std::streambuf(),
      pbuf_(nullptr),
      level_(log_level),
      stream_list_(),
      incomplete_line_(),
      log_cache_counter_(0),
      log_cache_(),
      log_time_cache_()
    {
      pbuf_ = new char[BUFFER_LENGTH];
      std::streambuf::setp(pbuf_, pbuf_ + BUFFER_LENGTH - 1);
    }

    LogStreamBuf::~LogStreamBuf()
    {
      sync();
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        clearCache();
        if (incomplete_line_.size() > 0)
          distribute_(incomplete_line_);
        delete[] pbuf_;
        pbuf_ = nullptr;
      }
    }

    int LogStreamBuf::overflow(int c)
    {
      if (c != traits_type::eof())
      {
        *pptr() = c;
        pbump(1);
        sync();
        return c;
      }
      else
      {
        return traits_type::eof();
      }
    }

    LogStreamBuf * LogStream::rdbuf()
    {
      return (LogStreamBuf *)std::ios::rdbuf();
    }

    LogStreamBuf * LogStream::operator->()
    {
      return rdbuf();
    }

    void LogStream::setLevel(std::string level)
    {
      if (rdbuf() == nullptr)
      {
        return;
      }

      // set the new level
      rdbuf()->level_ = level;
    }

    std::string LogStream::getLevel()
    {
      if (rdbuf() != nullptr)
      {
        return rdbuf()->level_;
      }
      else
      {
        return LogStreamBuf::UNKNOWN_LOG_LEVEL;
      }
    }

    // caching methods
    Size LogStreamBuf::getNextLogCounter_()
    {
      return ++log_cache_counter_;
    }

    bool LogStreamBuf::isInCache_(std::string const & line)
    {
      //cout << "LogCache (count)" << log_cache_.count(line) << endl;
      if (log_cache_.count(line) == 0)
      {
        return false;
      }
      else
      {
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

      if (log_cache_.size() > 1) // check if we need to remove one of the entries
      {
        // get smallest key
        map<Size, string>::iterator it = log_time_cache_.begin();

        // check if message occurred more then once
        if (log_cache_[it->second].counter != 0)
        {
          std::stringstream stream;
          stream << "<" << it->second << "> occurred " << ++log_cache_[it->second].counter << " times";
          extra_message = stream.str();
        }

        log_cache_.erase(it->second);
        log_time_cache_.erase(it);
      }

      Size counter_value = getNextLogCounter_();
      log_cache_[line].counter = 0;
      log_cache_[line].timestamp = counter_value;

      log_time_cache_[counter_value] = line;

      return extra_message;
    }

    void LogStreamBuf::clearCache()
    {
      // if there are any streams in our list, we
      // copy the line into that streams, too and flush them
      map<std::string, LogCacheStruct>::iterator it = log_cache_.begin();

      for (; it != log_cache_.end(); ++it)
      {
        if ((it->second).counter != 0)
        {
          std::stringstream stream;
          stream << "<" << it->first << "> occurred " << ++(it->second).counter << " times";
          distribute_(stream.str());
        }
      }
      // remove all entries from cache
      log_cache_.clear();
      log_time_cache_.clear();
    }

    int LogStreamBuf::sync()
    {
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        // sync our streambuffer...
        if (pptr() != pbase())
        {
          // check if we have attached streams, so we don't waste time to
          // prepare the output
          if (!stream_list_.empty())
          {
            char * line_start = pbase();
            char * line_end = pbase();

            static char buf[BUFFER_LENGTH];

            while (line_end < pptr())
            {
              // search for the first end of line
              for (; line_end < pptr() && *line_end != '\n'; line_end++)
              {
              }

              if (line_end >= pptr())
              {
                // Copy the incomplete line to the incomplete_line_ buffer
                size_t length = line_end - line_start;
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
                std::string outstring;
                std::swap(outstring, incomplete_line_); // init outstring, while resetting incomplete_line_ 
                outstring += &(buf[0]);

                // avoid adding empty lines to the cache
                if (outstring.empty())
                {
                  distribute_(outstring);
                }
                // check if we have already seen this log message
                else if (!isInCache_(outstring))
                {
                  // add line to the log cache
                  std::string extra_message = addToCache_(outstring);

                  // send outline (and extra_message) to attached streams
                  if (!extra_message.empty())
                    distribute_(extra_message);

                  distribute_(outstring);
                }

                // update the line pointers (increment both)
                line_start = ++line_end;
              }
            }
          }
          // remove all processed lines from the buffer
          pbump((int)(pbase() - pptr()));
        }

      } // ! OMP

      return 0;
    }

    void LogStreamBuf::distribute_(std::string outstring)
    {
      // if there are any streams in our list, we
      // copy the line into that streams, too and flush them
      std::list<StreamStruct>::iterator list_it = stream_list_.begin();
      for (; list_it != stream_list_.end(); ++list_it)
      {
        *(list_it->stream) << expandPrefix_(list_it->prefix, time(nullptr)).c_str()
                           << outstring.c_str() << std::endl;

        if (list_it->target != nullptr)
        {
          list_it->target->logNotify();
        }
      }
    }

    string LogStreamBuf::expandPrefix_
      (const std::string & prefix, time_t time) const
    {
      string::size_type   index = 0;
      Size copied_index = 0;
      string result("");

      while ((index = prefix.find("%", index)) != String::npos)
      {
        // append any constant parts of the string to the result
        if (copied_index < index)
        {
          result.append(prefix.substr(copied_index, index - copied_index));
          copied_index = (SignedSize)index;
        }

        if (index < prefix.size())
        {
          char    buffer[64] = "";
          char * buf = &(buffer[0]);

          switch (prefix[index + 1])
          {
          case '%':           // append a '%' (escape sequence)
            result.append("%");
            break;

          case 'y':           // append the message type (error/warning/information)
            result.append(level_);
            break;

          case 'T':           // time: HH:MM:SS
            strftime(buf, 64, "%H:%M:%S", localtime(&time));
            result.append(buf);
            break;

          case 't':           // time: HH:MM
            strftime(buf, 64, "%H:%M", localtime(&time));
            result.append(buf);
            break;

          case 'D':           // date: DD.MM.YYYY
            strftime(buf, 64, "%Y/%m/%d", localtime(&time));
            result.append(buf);
            break;

          case 'd':           // date: DD.MM.
            strftime(buf, 64, "%m/%d", localtime(&time));
            result.append(buf);
            break;

          case 'S':           // time+date: DD.MM.YYYY, HH:MM:SS
            strftime(buf, 64, "%Y/%m/%d, %H:%M:%S", localtime(&time));
            result.append(buf);
            break;

          case 's':           // time+date: DD.MM., HH:MM
            strftime(buf, 64, "%m/%d, %H:%M", localtime(&time));
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

    LogStreamNotifier::LogStreamNotifier() :
      registered_at_(nullptr)
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
      if (registered_at_ == nullptr)
        return;

      registered_at_->remove(stream_);
      registered_at_ = nullptr;
    }

    void LogStreamNotifier::registerAt(LogStream & log)
    {
      unregister();

      registered_at_ = &log;
      log.insertNotification(stream_, *this);
    }

    // keep the given buffer
    LogStream::LogStream(LogStreamBuf * buf, bool delete_buf, std::ostream * stream) :
      std::ios(buf),
      std::ostream(buf),
      delete_buffer_(delete_buf)
    {
      if (stream != nullptr)
      {
        insert(*stream);
      }
    }

    LogStream::~LogStream()
    {
      if (delete_buffer_)
      {
        // delete the streambuffer
        delete rdbuf();
        // set it to 0
        std::ios(nullptr);
      }
    }

    void LogStream::insert(std::ostream & stream)
    {
      if (!bound_() || hasStream_(stream))
      {
        return;
      }

      // we didn't find it - create a new entry in the list
      LogStreamBuf::StreamStruct s_struct;
      s_struct.stream = &stream;
      rdbuf()->stream_list_.push_back(s_struct);
    }

    void LogStream::remove(std::ostream & stream)
    {
      if (!bound_())
        return;

      StreamIterator it = findStream_(stream);
      if (it != rdbuf()->stream_list_.end())
      {
        rdbuf()->sync();
        // HINT: we do NOT clear the cache (because we cannot access it from here)
        //       and we do not flush incomplete_line_!!!
        rdbuf()->stream_list_.erase(it);
      }
    }

    void LogStream::insertNotification(std::ostream & s, LogStreamNotifier & target)
    {
      if (!bound_())
        return;

      insert(s);

      StreamIterator it = findStream_(s);
      (*it).target = &target;
    }

    LogStream::StreamIterator LogStream::findStream_(const std::ostream & s)
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

    bool LogStream::hasStream_(std::ostream & stream)
    {
      if (!bound_())
        return false;

      return findStream_(stream) != rdbuf()->stream_list_.end();
    }

    void LogStream::setPrefix(const std::ostream & s, const string & prefix)
    {
      if (!bound_())
        return;

      StreamIterator it = findStream_(s);
      if (it != rdbuf()->stream_list_.end())
      {
        (*it).prefix = prefix;
      }
    }

    void LogStream::setPrefix(const string & prefix)
    {
      if (!bound_())
        return;

      for (StreamIterator it = rdbuf()->stream_list_.begin(); it != rdbuf()->stream_list_.end(); ++it)
      {
        (*it).prefix = prefix;
      }
    }

    void LogStream::flush()
    {
      std::ostream::flush();
    }

    bool LogStream::bound_() const
    {
      LogStream * non_const_this = const_cast<LogStream *>(this);

      return non_const_this->rdbuf() != nullptr;
    }

  }   // namespace Logger

  // global StreamHandler
  OPENMS_DLLAPI StreamHandler STREAM_HANDLER;

  // global default logstream
  OPENMS_DLLAPI Logger::LogStream Log_fatal(new Logger::LogStreamBuf("FATAL_ERROR"), true, &cerr);
  OPENMS_DLLAPI Logger::LogStream Log_error(new Logger::LogStreamBuf("ERROR"), true, &cerr);
  OPENMS_DLLAPI Logger::LogStream Log_warn(new Logger::LogStreamBuf("WARNING"), true, &cout);
  OPENMS_DLLAPI Logger::LogStream Log_info(new Logger::LogStreamBuf("INFO"), true, &cout);
  // LOG_DEBUG is disabled by default, but will be enabled in TOPPAS.cpp or TOPPBase.cpp if started in debug mode (--debug or -debug X)
  OPENMS_DLLAPI Logger::LogStream Log_debug(new Logger::LogStreamBuf("DEBUG"), false); // last param should be 'true', but segfaults...

} // namespace OpenMS
