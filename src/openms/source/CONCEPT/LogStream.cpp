// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <ctime>
#include <mutex>
#include <vector>
#include <algorithm>    // std::min
#ifdef __GLIBCXX__
#include <cxxabi.h>     // abi::__forced_unwind
#endif
#include <OpenMS/CONCEPT/Colorizer.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/StreamHandler.h>

#include <sstream>
#include <iostream>

#define BUFFER_LENGTH 32768

using namespace std;

namespace
{
  /// Global mutex that serializes the final sink writes in
  /// OpenMS::Logger::LogStreamBuf::distribute_(). Multiple thread-local
  /// LogStreamBuf instances legitimately share the same destination ostream
  /// (e.g. std::cerr/std::cout) AND the same Colorizer of their level
  /// (see logColorizer_()), so both the stream writes and the Colorizer's
  /// internal state mutation must be serialized across threads (issue #9515).
  /// Intentionally heap-allocated and never freed so it stays valid even while
  /// the global LogStream objects are destroyed during static teardown; unlike
  /// the former OpenMP critical section, a plain std::mutex does not depend on
  /// the OpenMP runtime.
  std::mutex& logSinkMutex_()
  {
    static std::mutex* instance = new std::mutex();
    return *instance;
  }

  /// The Colorizer of the log levels in @p color, created on first use and never destroyed, like the global log
  /// streams. Not the public OpenMS::red, yellow or magenta: static initializers and destructors of other files may log
  /// before those are constructed or after they are destroyed. Used under logSinkMutex_() only.
  template<OpenMS::ConsoleColor color>
  OpenMS::Colorizer* logColorizer_()
  {
    static OpenMS::Colorizer* const instance = new OpenMS::Colorizer(color);
    return instance;
  }

  /// A new StreamStruct::id, unique for the process
  OpenMS::Size nextDestinationId_()
  {
    static std::atomic<OpenMS::Size> next_id{1};
    return next_id.fetch_add(1, std::memory_order_relaxed);
  }

  /// Thread-safe local-time conversion. std::localtime returns a pointer to a
  /// single process-wide static std::tm, which races across threads; the
  /// reentrant localtime_r/localtime_s write into a caller-provided struct.
  std::tm toLocalTime_(std::time_t t)
  {
    std::tm out{};
#ifdef OPENMS_WINDOWSPLATFORM
    localtime_s(&out, &t);
#else
    localtime_r(&t, &out);
#endif
    return out;
  }
}

namespace OpenMS
{
  namespace Logger
  {

    const time_t LogStreamBuf::MAX_TIME = numeric_limits<time_t>::max();
    const std::string LogStreamBuf::UNKNOWN_LOG_LEVEL = "UNKNOWN_LOG_LEVEL";

    LogStreamBuf::LogStreamBuf(const std::string& log_level, Colorizer* col)
      : std::streambuf(),
        level_(log_level),
        colorizer_(col)
    {
      pbuf_ = new char[BUFFER_LENGTH];
      std::streambuf::setp(pbuf_, pbuf_ + BUFFER_LENGTH - 1);
    }

    LogStreamBuf::LogStreamBuf(LogStreamBuf* source_buf, Colorizer* col)
      : std::streambuf(),
        level_(source_buf ? source_buf->level_ : UNKNOWN_LOG_LEVEL),
        colorizer_(col)
    {
      pbuf_ = new char[BUFFER_LENGTH];
      std::streambuf::setp(pbuf_, pbuf_ + BUFFER_LENGTH - 1);
      // Follow the source buffer's destinations, including later changes (see updateFromParent_())
      if (source_buf)
      {
        parent_ = source_buf;
        std::lock_guard<std::mutex> lock(logSinkMutex_());
        rebuildStreamList_();
      }
    }

    std::list<LogStreamBuf::StreamStruct>& LogStreamBuf::getStreamList_()
    {
      return stream_list_;
    }

    const std::list<LogStreamBuf::StreamStruct>& LogStreamBuf::getStreamList_() const
    {
      return stream_list_;
    }

    void LogStreamBuf::updateFromParent_()
    {
      if (parent_ != nullptr && parent_->version_.load(std::memory_order_acquire) != parent_version_)
      {
        std::lock_guard<std::mutex> lock(logSinkMutex_());
        updateFromParentLocked_();
      }
    }

    void LogStreamBuf::updateFromParentLocked_()
    {
      if (parent_ != nullptr && parent_->version_.load(std::memory_order_acquire) != parent_version_)
      {
        rebuildStreamList_();
      }
    }

    void LogStreamBuf::rebuildStreamList_()
    {
      // The sink mutex, held by the caller, also guards parent_->stream_list_ against concurrent changes.
      const std::vector<Size> previous = destinations_();
      parent_version_ = parent_->version_.load(std::memory_order_relaxed);
      stream_list_.clear();
      for (const StreamStruct& s : parent_->stream_list_)
      {
        const bool hidden = std::find(hidden_ids_.begin(), hidden_ids_.end(), s.id) != hidden_ids_.end();
        const bool own = std::any_of(own_streams_.begin(), own_streams_.end(),
                                     [&s](const StreamStruct& o) { return o.stream == s.stream; });
        if (!hidden && !own && !guarded_(s.stream))
        {
          stream_list_.push_back(s);
        }
      }
      for (const StreamStruct& s : own_streams_)
      {
        if (!guarded_(s.stream))
        {
          stream_list_.push_back(s);
        }
      }
      version_.fetch_add(1, std::memory_order_release);
      forgetRepeatsAfterChange_(previous);
    }

    bool LogStreamBuf::guarded_(const std::ostream* stream) const
    {
      return std::find(guarded_streams_.begin(), guarded_streams_.end(), stream) != guarded_streams_.end();
    }

    std::list<LogStreamBuf::StreamStruct>& LogStreamBuf::prefixableStreams_()
    {
      return parent_ == nullptr ? stream_list_ : own_streams_;
    }

    void LogStreamBuf::destinationsChanged_()
    {
      if (parent_ != nullptr)
      {
        rebuildStreamList_();
      }
      else
      {
        version_.fetch_add(1, std::memory_order_release);
      }
    }

    std::vector<Size> LogStreamBuf::destinations_() const
    {
      std::vector<Size> ids;
      for (const StreamStruct& s : stream_list_)
      {
        ids.push_back(s.id);
      }
      std::sort(ids.begin(), ids.end());
      return ids;
    }

    void LogStreamBuf::forgetRepeatsAfterChange_(const std::vector<Size>& previous)
    {
      if (log_cache_.empty() || previous == destinations_())
      {
        return;
      }
      // The cache stands for the messages that the previous destinations got. Their pending repeat counts go to
      // those that remain (a removed one may be destroyed already); a new one gets every message from now on.
      for (const std::string& summary : repeatSummaries_())
      {
        for (StreamStruct& s : stream_list_)
        {
          if (std::binary_search(previous.begin(), previous.end(), s.id))
          {
            write_(s, summary);
          }
        }
      }
      log_cache_.clear();
      log_time_cache_.clear();
    }

    std::vector<std::string> LogStreamBuf::repeatSummaries_() const
    {
      std::vector<std::string> summaries;
      for (const auto& [line, entry] : log_cache_)
      {
        if (entry.counter != 0)
        {
          summaries.push_back("<" + line + "> occurred " + std::to_string(entry.counter + 1) + " times");
        }
      }
      return summaries;
    }

    LogStreamBuf::~LogStreamBuf()
    {
      // Flush whatever is left. distribute_() serializes on a global std::mutex
      // that is intentionally leaked (never destroyed), so it stays valid while
      // the global LogStream objects are torn down at static destruction. Unlike
      // the former OpenMP critical section, a std::mutex does not depend on the
      // OpenMP runtime, so locking here during teardown is safe (issue #9515).
      syncLF_();
      {
        clearCache();
        if (!incomplete_line_.empty())
        {
          distribute_(incomplete_line_);
        }
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
      rdbuf()->level_ = std::move(level);
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
      if (!log_cache_.contains(line))
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
      std::string extra_message;
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
      std::lock_guard<std::mutex> lock(logSinkMutex_());
      // Apply changes of the parent's destinations first: the repeat counts pending from before go to the
      // destinations that got the repeated messages (see forgetRepeatsAfterChange_()).
      updateFromParentLocked_();
      for (const std::string& summary : repeatSummaries_())
      {
        for (StreamStruct& s : stream_list_)
        {
          write_(s, summary);
        }
      }
      log_cache_.clear();
      log_time_cache_.clear();
    }

    void LogStreamBuf::distribute_(const std::string& outstring)
    {
      // Serialize the final writes across threads. Multiple thread-local
      // LogStreamBuf instances legitimately share the same destination ostream
      // (e.g. std::cerr/std::cout) AND the same Colorizer of their level,
      // so both the stream writes and the Colorizer's
      // internal state mutation must be serialized (issue #9515).
      std::lock_guard<std::mutex> lock(logSinkMutex_());

      // Pick up changes of the parent's destinations while holding the lock, so that a stream removed
      // there (and possibly destroyed afterwards, e.g. by StreamHandler) is never written to.
      updateFromParentLocked_();

      for (StreamStruct& s : stream_list_)
      {
        write_(s, outstring);
      }
    }

    void LogStreamBuf::distributeNew_(const std::string& line)
    {
      std::lock_guard<std::mutex> lock(logSinkMutex_());
      // Apply changes of the parent's destinations before evicting. If the destinations changed, the pending repeat
      // counts went to the destinations that got the repeated messages, and the cache is empty now (see
      // forgetRepeatsAfterChange_()). So the count of a message evicted here belongs to the current destinations:
      // a destination inserted since the message was logged never gets it.
      updateFromParentLocked_();
      const std::string extra_message = addToCache_(line);
      if (!extra_message.empty())
      {
        for (StreamStruct& s : stream_list_)
        {
          write_(s, extra_message);
        }
      }
      for (StreamStruct& s : stream_list_)
      {
        write_(s, line);
      }
    }

    void LogStreamBuf::write_(StreamStruct& s, const std::string& line)
    {
      try
      {
        if (colorizer_)
        {
          *(s.stream) << (*colorizer_)(); // enable color
        }

        *(s.stream) << expandPrefix_(s.prefix, time(nullptr)) << line;

        if (colorizer_)
        {
          *(s.stream) << (*colorizer_).undo(); // disable color
        }
        *(s.stream) << std::endl;
      }
#ifdef __GLIBCXX__
      catch (abi::__forced_unwind&)
      {
        throw; // thread cancellation, as std::ostream itself does
      }
#endif
      catch (...)
      {
        // A destination with exceptions enabled failed (e.g. std::ios_base::failure). It must not keep the line from
        // the other destinations, nor make logging (or a change of the destinations, which writes repeat counts) throw.
      }
    }

    int LogStreamBuf::syncLF_()
    {
      // sync our stream buffer...
      if (pptr() != pbase())
      {
        updateFromParent_();
        // check if we have attached streams, so we don't waste time to
        // prepare the output
        if (!stream_list_.empty())
        {
          char *line_start = pbase();
          char *line_end = pbase();

          while (line_end < pptr())
          {
            // search for the first end of line
            for (; line_end < pptr() && *line_end != '\n'; line_end++)
            {
            }

            if (line_end >= pptr())
            {
              // No newline yet: stash the partial line for the next flush.
              // Appending straight into incomplete_line_ (a std::string) avoids
              // the former shared function-local 'static' scratch buffer, which
              // is raced on by concurrent thread-local LogStreamBufs (issue #9515).
              incomplete_line_.append(line_start, (size_t) (line_end - line_start));

              // mark everything as read
              line_end = pptr() + 1;
            }
            else
            {
              // A full line ends at line_end (the '\n'); [line_start, line_end) is
              // its content without the newline. Assemble the string to be written,
              // prepending any leftover from the previous flush (incomplete_line_).
              std::string outstring;
              std::swap(outstring, incomplete_line_); // init outstring, while resetting incomplete_line_
              outstring.append(line_start, (size_t) (line_end - line_start));

              // avoid adding empty lines to the cache
              if (outstring.empty())
              {
                distribute_(outstring);
              }
                // check if we have already seen this log message
              else if (!isInCache_(outstring))
              {
                // add line to the log cache and send it (after the repeat count of an evicted message) to attached streams
                distributeNew_(outstring);
              }

              // update the line pointers (increment both)
              line_start = ++line_end;
            }
          }
        }
        // remove all processed lines from the buffer
        pbump((int) (pbase() - pptr()));
      }
      return 0;
    }

    int LogStreamBuf::sync()
    {
      int ret = 0;

        ret = syncLF_();
      
      return ret;
    }

    string LogStreamBuf::expandPrefix_
      (const std::string & prefix, time_t time) const
    {
      string::size_type   index = 0;
      Size copied_index = 0;
      string result;

      while ((index = prefix.find('%', index)) != std::string::npos)
      {
        // append any constant parts of the string to the result
        if (copied_index < index)
        {
          result.append(StringUtils::substr(prefix, copied_index, index - copied_index));
          copied_index = (SignedSize)index;
        }

        if (index < prefix.size())
        {
          char    buffer[64] = "";
          char * buf = &(buffer[0]);
          // reentrant local-time; std::localtime would share a static std::tm
          [[maybe_unused]] std::tm tmv = toLocalTime_(time);

          switch (prefix[index + 1])
          {
          case '%':           // append a '%' (escape sequence)
            result.append("%");
            break;

          case 'y':           // append the message type (error/warning/information)
            result.append(level_);
            break;

          case 'T':           // time: HH:MM:SS
            strftime(buf, 64, "%H:%M:%S", &tmv);
            result.append(buf);
            break;

          case 't':           // time: HH:MM
            strftime(buf, 64, "%H:%M", &tmv);
            result.append(buf);
            break;

          case 'D':           // date: DD.MM.YYYY
            strftime(buf, 64, "%Y/%m/%d", &tmv);
            result.append(buf);
            break;

          case 'd':           // date: DD.MM.
            strftime(buf, 64, "%m/%d", &tmv);
            result.append(buf);
            break;

          case 'S':           // time+date: DD.MM.YYYY, HH:MM:SS
            strftime(buf, 64, "%Y/%m/%d, %H:%M:%S", &tmv);
            result.append(buf);
            break;

          case 's':           // time+date: DD.MM., HH:MM
            strftime(buf, 64, "%m/%d, %H:%M", &tmv);
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
        result.append(StringUtils::substr(prefix, copied_index, prefix.size() - copied_index));
      }

      return result;
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
        // delete the stream buffer
        delete rdbuf();
        std::ios::rdbuf(nullptr);
      }
    }

    // Changes of the destinations hold the sink mutex: following (thread-local) buffers read their
    // parent's list under it, and distribute_() writes to the destinations under it.

    void LogStream::insert(std::ostream & stream)
    {
      if (!bound_())
      {
        return;
      }
      followParent_();
      if (hasStream_(stream))
      {
        return;
      }
      LogStreamBuf* buf = rdbuf();
      // lines and repeat counts pending from before belong to the current destinations, not to the new one
      buf->sync();
      buf->clearCache();
      std::lock_guard<std::mutex> lock(logSinkMutex_());
      bool inherited = false;
      if (buf->parent_ != nullptr)
      {
        // undo the removal of a destination that the parent still has
        for (const LogStreamBuf::StreamStruct& s : buf->parent_->stream_list_)
        {
          if (s.stream == &stream)
          {
            buf->hidden_ids_.erase(std::remove(buf->hidden_ids_.begin(), buf->hidden_ids_.end(), s.id), buf->hidden_ids_.end());
            inherited = true;
          }
        }
      }
      std::list<LogStreamBuf::StreamStruct>& entries = buf->parent_ == nullptr ? buf->stream_list_ : buf->own_streams_;
      // an own destination may be there already, suppressed by a LogSinkGuard
      const bool present = std::any_of(entries.begin(), entries.end(),
                                       [&stream](const LogStreamBuf::StreamStruct& s) { return s.stream == &stream; });
      if (!inherited && !present)
      {
        LogStreamBuf::StreamStruct s_struct;
        s_struct.stream = &stream;
        s_struct.id = nextDestinationId_();
        entries.push_back(s_struct);
      }
      buf->destinationsChanged_();
      updateState_();
    }

    void LogStream::remove(std::ostream & stream)
    {
      if (!bound_())
        return;

      followParent_();
      LogStreamBuf* buf = rdbuf();
      if (!hasStream_(stream) && !buf->guarded_(&stream)) // a removal in a guarded scope lasts beyond it
      {
        return;
      }
      buf->sync();
      // pending repeat counts also go to the stream that is removed; an incomplete line stays for the remaining ones
      buf->clearCache();
      std::lock_guard<std::mutex> lock(logSinkMutex_());
      detachLocked_(stream);
      buf->destinationsChanged_();
      updateState_();
    }

    LogStream::GuardedRemoval LogStream::detachLocked_(std::ostream& stream)
    {
      GuardedRemoval removal;
      LogStreamBuf* buf = rdbuf();
      std::list<LogStreamBuf::StreamStruct>& entries = buf->parent_ == nullptr ? buf->stream_list_ : buf->own_streams_;
      for (const LogStreamBuf::StreamStruct& s : entries)
      {
        if (s.stream == &stream)
        {
          removal.reinsert = true;
          removal.prefix = s.prefix;
          removal.id = s.id;
        }
      }
      entries.remove_if([&stream](const LogStreamBuf::StreamStruct& s) { return s.stream == &stream; });
      if (buf->parent_ != nullptr)
      {
        // hide the parent's destination from this buffer
        for (const LogStreamBuf::StreamStruct& s : buf->parent_->stream_list_)
        {
          if (s.stream == &stream && std::find(buf->hidden_ids_.begin(), buf->hidden_ids_.end(), s.id) == buf->hidden_ids_.end())
          {
            buf->hidden_ids_.push_back(s.id);
          }
        }
      }
      return removal;
    }

    LogStream::GuardedRemoval LogStream::removeForGuard_(std::ostream& stream)
    {
      if (!bound_())
      {
        return {};
      }
      followParent_();
      drain_();
      if (LogStream* follower = threadLocalFollower_())
      {
        follower->drain_();
      }
      LogStreamBuf* buf = rdbuf();
      std::lock_guard<std::mutex> lock(logSinkMutex_());
      buf->updateFromParentLocked_();
      GuardedRemoval removal;
      if (buf->parent_ == nullptr)
      {
        removal = detachLocked_(stream);
      }
      else
      {
        // suppressed on this thread until the guard ends, also if the global stream inserts it again meanwhile
        buf->guarded_streams_.push_back(&stream);
      }
      buf->destinationsChanged_();
      updateState_();
      return removal;
    }

    void LogStream::restoreForGuard_(std::ostream& stream, const GuardedRemoval& removal)
    {
      if (!bound_())
      {
        return;
      }
      followParent_();
      // text written in the guarded scope goes to the destinations without the stream
      drain_();
      if (LogStream* follower = threadLocalFollower_())
      {
        follower->drain_();
      }
      LogStreamBuf* buf = rdbuf();
      std::lock_guard<std::mutex> lock(logSinkMutex_());
      buf->updateFromParentLocked_();
      if (buf->parent_ != nullptr)
      {
        auto guarded = std::find(buf->guarded_streams_.begin(), buf->guarded_streams_.end(), &stream);
        if (guarded != buf->guarded_streams_.end())
        {
          buf->guarded_streams_.erase(guarded); // one entry: an enclosing guard on the same stream keeps it suppressed
        }
      }
      std::list<LogStreamBuf::StreamStruct>& entries = buf->stream_list_;
      const bool present = std::any_of(entries.begin(), entries.end(),
                                       [&stream](const LogStreamBuf::StreamStruct& s) { return s.stream == &stream; });
      if (buf->parent_ == nullptr && removal.reinsert && !present)
      {
        LogStreamBuf::StreamStruct s_struct;
        s_struct.stream = &stream;
        s_struct.prefix = removal.prefix;
        // the same destination as before, so that threads that hid it on their thread-local stream keep it hidden
        s_struct.id = removal.id;
        entries.push_back(s_struct);
      }
      buf->destinationsChanged_();
      updateState_();
    }

    void LogStream::drain_()
    {
      flushIncomplete();
      rdbuf()->clearCache();
    }

    LogStream* LogStream::threadLocalFollower_()
    {
      if (rdbuf()->parent_ != nullptr)
      {
        return nullptr;
      }
      if (this == &getGlobalLogFatal()) return &getThreadLocalLogFatal();
      if (this == &getGlobalLogError()) return &getThreadLocalLogError();
      if (this == &getGlobalLogWarn()) return &getThreadLocalLogWarn();
      if (this == &getGlobalLogInfo()) return &getThreadLocalLogInfo();
      if (this == &getGlobalLogDebug()) return &getThreadLocalLogDebug();
      return nullptr;
    }

    void LogStream::removeAllStreams()
    {
      if (!bound_())
        return;

      LogStreamBuf* buf = rdbuf();
      buf->sync();
      // Distribute any incomplete line before clearing streams
      if (!buf->incomplete_line_.empty())
      {
        buf->distribute_(buf->incomplete_line_);
        buf->incomplete_line_.clear();
      }
      buf->clearCache();
      std::lock_guard<std::mutex> lock(logSinkMutex_());
      // Flush all streams before clearing the list. Under the lock, like the writes in distribute_(),
      // and after updating from the parent, which may have removed (and destroyed) a stream.
      buf->updateFromParentLocked_();
      for (auto& stream_struct : buf->stream_list_)
      {
        if (stream_struct.stream != nullptr)
        {
          stream_struct.stream->flush();
        }
      }
      buf->stream_list_.clear();
      if (buf->parent_ != nullptr)
      {
        // hide the parent's current destinations from this buffer
        buf->own_streams_.clear();
        for (const LogStreamBuf::StreamStruct& s : buf->parent_->stream_list_)
        {
          if (std::find(buf->hidden_ids_.begin(), buf->hidden_ids_.end(), s.id) == buf->hidden_ids_.end())
          {
            buf->hidden_ids_.push_back(s.id);
          }
        }
      }
      buf->destinationsChanged_();
      updateState_();
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
      {
        return false;
      }
      return findStream_(stream) != rdbuf()->stream_list_.end();
    }

    void LogStream::setPrefix(const std::ostream & s, const string & prefix)
    {
      if (!bound_())
      {
        return;
      }
      LogStreamBuf* buf = rdbuf();
      std::lock_guard<std::mutex> lock(logSinkMutex_());
      buf->updateFromParentLocked_(); // pending repeat counts keep the prefix of the lines they count
      std::list<LogStreamBuf::StreamStruct>& entries = buf->prefixableStreams_();
      auto entry = std::find_if(entries.begin(), entries.end(),
                                [&s](const LogStreamBuf::StreamStruct& e) { return e.stream == &s; });
      if (entry != entries.end())
      {
        entry->prefix = prefix;
        buf->destinationsChanged_();
      }
    }

    void LogStream::setPrefix(const string & prefix)
    {
      if (!bound_())
      {
        return;
      }
      LogStreamBuf* buf = rdbuf();
      std::lock_guard<std::mutex> lock(logSinkMutex_());
      buf->updateFromParentLocked_(); // pending repeat counts keep the prefix of the lines they count
      for (LogStreamBuf::StreamStruct& entry : buf->prefixableStreams_())
      {
        entry.prefix = prefix;
      }
      buf->destinationsChanged_();
    }

    bool LogStream::bound_() const
    {
      LogStream * non_const_this = const_cast<LogStream *>(this);

      return non_const_this->rdbuf() != nullptr;
    }

    void LogStream::followParent_()
    {
      if (!bound_())
      {
        return;
      }
      rdbuf()->updateFromParent_();
      updateState_();
    }

    void LogStream::updateState_()
    {
      const LogStreamBuf* buf = rdbuf();
      if (buf->parent_ == nullptr)
      {
        return;
      }
      if (buf->stream_list_.empty())
      {
        setstate(std::ios_base::badbit);
      }
      else if (rdstate() != std::ios_base::goodbit)
      {
        clear();
      }
    }

    /// Lets the thread-local accessors update their stream on each use
    struct ThreadLocalLogAccess
    {
      static LogStream& use(LogStream& stream)
      {
        stream.followParent_();
        return stream;
      }
    };

    void LogStream::flush()
    {
      std::ostream::flush();
    }

    bool LogStream::hasStream(std::ostream & stream)
    {
      if (!bound_())
      {
        return false;
      }
      rdbuf()->updateFromParent_();
      return hasStream_(stream);
    }

    void LogStream::flushIncomplete()
    {
      if (!bound_())
        return;

      rdbuf()->sync();
      // Distribute any incomplete line (text not terminated by newline)
      if (!rdbuf()->incomplete_line_.empty())
      {
        rdbuf()->distribute_(rdbuf()->incomplete_line_);
        rdbuf()->incomplete_line_.clear();
      }
    }

  }   // namespace Logger


  // global StreamHandler
  OPENMS_DLLAPI StreamHandler STREAM_HANDLER;

  //
  // Global log stream accessor functions (for configuration purposes)
  // WARNING: Direct logging to these streams is NOT thread-safe.
  // Use OPENMS_LOG_* macros (which use thread-local streams) for actual logging.
  //
  // Each global stream is created on its first use, so that static initializers of other files can log, and is never
  // destroyed: thread-local streams follow them and may still log, or flush when they are destroyed, after static
  // destruction has begun. Their own pending output is flushed at exit instead (see GlobalLogStreamsFlusher).
  //
  Logger::LogStream& getGlobalLogFatal()
  {
    static Logger::LogStream& stream = *new Logger::LogStream(new Logger::LogStreamBuf("FATAL_ERROR", logColorizer_<ConsoleColor::RED>()), true, &cerr);
    return stream;
  }

  Logger::LogStream& getGlobalLogError()
  {
    static Logger::LogStream& stream = *new Logger::LogStream(new Logger::LogStreamBuf("ERROR", logColorizer_<ConsoleColor::RED>()), true, &cerr);
    return stream;
  }

  Logger::LogStream& getGlobalLogWarn()
  {
    static Logger::LogStream& stream = *new Logger::LogStream(new Logger::LogStreamBuf("WARNING", logColorizer_<ConsoleColor::YELLOW>()), true, &cerr);
    return stream;
  }

  Logger::LogStream& getGlobalLogInfo()
  {
    static Logger::LogStream& stream = *new Logger::LogStream(new Logger::LogStreamBuf("INFO", nullptr), true, &cout);
    return stream;
  }

  Logger::LogStream& getGlobalLogDebug()
  {
    // OPENMS_LOG_DEBUG is disabled by default, but will be enabled in TOPPAS.cpp or TOPPBase.cpp if started in debug mode (--debug or -debug X)
    static Logger::LogStream& stream = *new Logger::LogStream(new Logger::LogStreamBuf("DEBUG", logColorizer_<ConsoleColor::MAGENTA>()), true);
    return stream;
  }

  namespace
  {
    /// Flushes the global streams at exit, as their destructors would (see above)
    struct GlobalLogStreamsFlusher
    {
      ~GlobalLogStreamsFlusher()
      {
        for (Logger::LogStream* log : {&getGlobalLogFatal(), &getGlobalLogError(), &getGlobalLogWarn(), &getGlobalLogInfo(), &getGlobalLogDebug()})
        {
          log->flushIncomplete();
          log->rdbuf()->clearCache();
        }
      }
    } global_log_streams_flusher;
  }

  //
  // Thread-local log stream accessors
  // Each thread gets its own LogStream instance with a private buffer that follows the output
  // destinations of the global instance, including later changes (see LogStreamBuf::updateFromParent_()).
  // Reconfigure the global stream to redirect or suppress output of all threads, the thread-local
  // stream returned here for the calling thread only. Each call applies changes of the global destinations,
  // so that a stream without destinations is in a failed state and does not format messages (see
  // LogStream::updateState_()).
  //
  namespace
  {
    /**
      The calling thread's stream that follows @p global. A thread destroys its thread_local objects when it ends; the
      main thread does so before static destructors and atexit handlers run. A message logged after that, e.g. from such
      a destructor or from the destructor of another thread_local object, goes to a stream that is created then and never
      destroyed, instead of the destroyed one. Lines written to it reach the destinations; repeat counts and an incomplete
      line pending at the end are not written.
    */
    template<int level>
    Logger::LogStream& threadLocalStream_(Logger::LogStream& global, Colorizer* color)
    {
      thread_local bool destroyed = false; // trivially destructible, so it stays valid while the thread ends
      struct Owner
      {
        Logger::LogStream stream;
        Owner(Logger::LogStream& g, Colorizer* c) : stream(new Logger::LogStreamBuf(g.rdbuf(), c), true) {}
        ~Owner() { destroyed = true; }
      };
      if (destroyed)
      {
        thread_local Logger::LogStream* late = nullptr;
        if (late == nullptr)
        {
          late = new Logger::LogStream(new Logger::LogStreamBuf(global.rdbuf(), color), true);
        }
        return Logger::ThreadLocalLogAccess::use(*late);
      }
      thread_local Owner owner(global, color);
      return Logger::ThreadLocalLogAccess::use(owner.stream);
    }
  }

  Logger::LogStream& getThreadLocalLogFatal() { return threadLocalStream_<0>(getGlobalLogFatal(), logColorizer_<ConsoleColor::RED>()); }
  Logger::LogStream& getThreadLocalLogError() { return threadLocalStream_<1>(getGlobalLogError(), logColorizer_<ConsoleColor::RED>()); }
  Logger::LogStream& getThreadLocalLogWarn() { return threadLocalStream_<2>(getGlobalLogWarn(), logColorizer_<ConsoleColor::YELLOW>()); }
  Logger::LogStream& getThreadLocalLogInfo() { return threadLocalStream_<3>(getGlobalLogInfo(), nullptr); }
  Logger::LogStream& getThreadLocalLogDebug() { return threadLocalStream_<4>(getGlobalLogDebug(), logColorizer_<ConsoleColor::MAGENTA>()); }

  void setConsoleDebugLogging(bool enabled)
  {
    Logger::LogStream& local_debug = getThreadLocalLogDebug();
    // Flush to the previous destinations so pending or cached messages do not cross the switch.
    local_debug.flushIncomplete();
    local_debug.rdbuf()->clearCache();
    if (enabled)
    {
      getGlobalLogDebug().insert(cout);
      local_debug.insert(cout);
    }
    else
    {
      getGlobalLogDebug().remove(cout);
      local_debug.remove(cout);
    }
  }

} // namespace OpenMS
