// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Stephan Aiche, Andreas Bertsch$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <atomic>
#include <sstream>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include <map>

namespace OpenMS
{
  class Colorizer;

  /**
    @brief Log streams

    Logging, filtering, and storing messages.
    Many programs emit warning messages, error messages, or simply
    information and remarks to their users. The  LogStream
    class provides a convenient and straight-forward interface
    to classify these messages according to their importance
    (via the loglevel), filter and store them in files or
    write them to streams. \par
    As the LogStream class is derived from ostream, it behaves
    as any ostream object. Additionally you may associate
    streams with each LogStream object that catch only
    messages of certain loglevels. So the user might decide to
    redirect all error messages to cerr, all warning messages
    to cout and all information to a file. \par
    Along with each message its time of creation and its loglevel
    is stored. So the user might also decide to store all
    errors he got in the last two hours or alike. \par
    The LogStream class heavily relies on the LogStreamBuf
    class, which does the actual buffering and storing, but is only
    of interest if you want to implement a derived class, as the
    actual user interface is implemented in the LogStream class.

    @ingroup Concept
  */
  namespace Logger
  {
    // forward declarations
    class LogStream;

    /**
      @brief Stream buffer used by LogStream.

      This class implements the low level behavior of
      LogStream . It takes care of the buffers and stores
      the lines written into the LogStream object.
      It also contains a list of streams that are associated with
      the LogStream object. This list contains pointers to the
      streams and their minimum and maximum log level.
      Each line entered in the LogStream is marked with its
      time (in fact, the time LogStreamBuf::sync was called) and its
      loglevel. The loglevel is determined by either the current
      loglevel (as set by  LogStream::setLevel or a temporary
      level (as set by LogStream::level for a single line only).
      For each line stored, the list of associated streams is checked
      whether the loglevel falls into the range declared by the
      stream's minimum and maximum level. If this condition is met,
      the logline (with its prefix, see  LogStream::setPrefix )
      is also copied to the associated stream and this stream is
      flushed, too.
    */
    class OPENMS_DLLAPI LogStreamBuf :
      public std::streambuf
    {

      friend class LogStream;

public:

      /// @name Constants
      //@{
      static const time_t MAX_TIME;
      static const std::string UNKNOWN_LOG_LEVEL;
      //@}

      /// @name Constructors and Destructors
      //@{


      /**
        Create a new LogStreamBuf object and set the level to @p log_level

        @param[in] log_level The log level of the LogStreamBuf (default is unknown)
        @param[in] col If messages should be colored, provide a colorizer here
      */
      LogStreamBuf(const std::string& log_level = UNKNOWN_LOG_LEVEL, Colorizer* col = nullptr);

      /**
        Create a LogStreamBuf that writes to the destinations of another LogStreamBuf.

        This constructor is used for thread-local logging: each thread has its own buffer (for
        thread safety) that follows the destinations of the global LogStream instance, including
        changes made after this buffer was created. Destinations inserted into or removed from the
        LogStream of this buffer only affect this buffer (see LogStream::insert()).

        @param[in] source_buf The LogStreamBuf whose destinations should be followed
        @param[in] col If messages should be colored, provide a colorizer here
      */
      LogStreamBuf(LogStreamBuf* source_buf, Colorizer* col = nullptr);

      /**
        Destruct the buffer and free all stored messages strings.
      */
      ~LogStreamBuf() override;

      //@}

      /// @name Stream methods
      //@{

      /**
        This method is called as soon as the ostream is flushed
        (especially this method is called by flush or endl).
        It transfers the contents of the streambufs putbuffer
        into a logline if a newline or linefeed character
        is found in the buffer ("\n" or "\r" resp.).
        The line is then removed from the putbuffer.
        Incomplete lines (not terminated by "\n" / "\r" are
        stored in incomplete_line_.
      */
      int sync() override;

      /**
        This method calls sync and <tt>streambuf::overflow(c)</tt> to
        prevent a buffer overflow.
      */
      int overflow(int c = -1) override;
      //@}


      /// @name Level methods
      //@{
      /**
        Set the level of the LogStream

        @param[in] level The new LogLevel
      */
      void setLevel(std::string level);


      /**
        Returns the LogLevel of this LogStream
      */
      std::string getLevel();
      //@}

      /**
        @brief Holds a stream that is connected to the LogStream.
        It also includes the minimum and maximum level at which the
        LogStream redirects messages to this stream.
      */
      struct OPENMS_DLLAPI StreamStruct
      {
        std::ostream * stream;
        std::string         prefix;
        /// Unique for each insertion: a stream inserted again, or a new stream at the address of a destroyed one, is a new destination
        Size id;

        StreamStruct() :
          stream(nullptr),
          id(0)
        {}
      };

      /**
        Checks if some of the cached entries where sent more then once
        to the LogStream and (if necessary) prints a corresponding messages
        into all affected Logs
      */
      void clearCache();

protected:

      /// Distribute a new message to connected streams.
      void distribute_(const std::string& outstring);

      /**
        Adds @p line, which is not in the cache, to the cache and distributes it, after the repeat count of a message it
        evicts. Under the sink mutex and after applying changes of parent_'s destinations, so that the count goes only
        to destinations that got the repeated message.
      */
      void distributeNew_(const std::string& line);

      /// Interpret the prefix format string and return the expanded prefix.
      std::string expandPrefix_(const std::string & prefix, time_t time) const;

      /// Returns the stream list (owned or shared)
      std::list<StreamStruct>& getStreamList_();
      const std::list<StreamStruct>& getStreamList_() const;

      /// Updates stream_list_ if parent_ changed its destinations since the last update. Locks the sink mutex.
      void updateFromParent_();

      /// Same as updateFromParent_(), but the caller holds the sink mutex.
      void updateFromParentLocked_();

      /// Rebuilds stream_list_ from parent_'s destinations and this buffer's own changes. The caller holds the sink mutex.
      void rebuildStreamList_();

      /// Is @p stream suppressed by a LogSinkGuard on this following buffer (see guarded_streams_)?
      bool guarded_(const std::ostream* stream) const;

      /// Entries whose prefix may be changed: stream_list_, or own_streams_ for a following buffer (the prefix of parent_'s destinations is set on parent_)
      std::list<StreamStruct>& prefixableStreams_();

      /// Publishes a change of the destinations to following buffers (or rebuilds stream_list_). The caller holds the sink mutex.
      void destinationsChanged_();

      char * pbuf_ = nullptr;
      std::string             level_;
      std::list<StreamStruct> stream_list_;  ///< Destinations of this buffer (for a following buffer: derived from parent_, own_streams_ and hidden_ids_)
      std::string             incomplete_line_;
      /// Buffer whose destinations this buffer follows (the global buffer of a thread-local one), or nullptr
      LogStreamBuf* parent_ = nullptr;
      /// Incremented whenever stream_list_ changes, so that following buffers know when to update
      std::atomic<Size> version_{0};
      /// version_ of parent_ that stream_list_ reflects
      Size parent_version_ = 0;
      /// Destinations inserted into this following buffer; they replace parent_'s entry for the same stream
      std::list<StreamStruct> own_streams_;
      /// Destinations of parent_ (their StreamStruct::id) removed from this following buffer. Kept when parent_ removes them:
      /// a LogSinkGuard on parent_ restores a destination with its id, and ids are never reused.
      std::vector<Size> hidden_ids_;
      /// Streams that a LogSinkGuard suppresses on this following buffer (one entry per guard), however they are inserted
      std::vector<const std::ostream*> guarded_streams_;
      Colorizer* colorizer_ = nullptr; ///< optional Colorizer to color the output to stdout/stdcerr (if attached)
      /// @name Caching
      //@{

      /**
        @brief Holds a counter of occurrences and an index for the occurrence sequence of the corresponding log message
      */
      struct LogCacheStruct
      {
        Size timestamp;
        int counter;
      };

      /**
        Sequential counter to remember the sequence of occurrence
        of the cached log messages
      */
      Size log_cache_counter_ = 0;

      /// Cache of the last two log messages
      std::map<std::string, LogCacheStruct> log_cache_;
      /// Cache of the occurrence sequence of the last two log messages
      std::map<Size, std::string> log_time_cache_;

      /// Checks if the line is already in the cache
      bool isInCache_(std::string const & line);

      /**
        Adds the new line to the cache and removes an old one
        if necessary

        @param[in] line The Log message that should be added to the cache
        @return An additional massage if a re-occurring message was removed
        from the cache
      */
      std::string addToCache_(std::string const & line);

      /// Returns the next free index for a log message
      Size getNextLogCounter_();

      /// Non-lock acquiring sync function called in the d'tor
      int syncLF_();

      /// "<message> occurred N times" for each cached message that was repeated
      std::vector<std::string> repeatSummaries_() const;

      /**
        After stream_list_ changed from the destinations @p previous (see destinations_()): writes the pending repeat
        counts to the destinations that remain and clears the cache, so that a new destination gets every message.
        The caller holds the sink mutex.
      */
      void forgetRepeatsAfterChange_(const std::vector<Size>& previous);
      //@}

      /// The destinations in stream_list_ (their StreamStruct::id), sorted
      std::vector<Size> destinations_() const;

      /// Writes @p line to the destination @p s, with its prefix and color. The caller holds the sink mutex.
      void write_(StreamStruct& s, const std::string& line);
    };

    /**
      @brief Log Stream Class.

      Defines a log stream which features a cache and some formatting.
      For the developer, however, only some macros are of interest which
      will push the message that follows them into the
      appropriate stream:

      Macros:
        - OPENMS_LOG_FATAL_ERROR
        - OPENMS_LOG_ERROR (non-fatal error are reported (processing continues))
        - OPENMS_LOG_WARN  (warning, a piece of information which should be read by the user, should be logged)
        - OPENMS_LOG_INFO (information, e.g. a status should be reported)
        - OPENMS_LOG_DEBUG (general debugging information -  output be written to cout if debug_level > 0)

      To use a specific logger of a log level simply use it as cerr or cout: <br>
      <code> OPENMS_LOG_ERROR << " A bad error occurred ..."  </code>
      <br>
      Which produces an error message in the log.

      @note The OPENMS_LOG_* macros are thread-safe: each thread logs through its
      own thread-local LogStream/LogStreamBuf (see getThreadLocalLog*()), so the
      per-stream buffers and caches are never shared. The final writes to the
      shared sink(s) (e.g. @c std::cerr / @c std::cout) and to the shared Colorizer
      are serialized by a global mutex inside LogStreamBuf::distribute_(). The global
      LogStream objects returned by getGlobalLog*() are NOT thread-safe for direct
      logging. Configure them from one thread at a time; every thread-local stream
      follows their destinations from its next message on.

    */
    class OPENMS_DLLAPI LogStream :
      public std::ostream
    {
public:

      /// @name Constructors and Destructors
      //@{

      /**
        Creates a new LogStream object that is not associated with any stream.
        If the argument <tt>stream</tt> is set to an output stream (e.g. <tt>cout</tt>)
        all output is send to that stream.

        @param[in]	buf
        @param[in]  delete_buf
        @param[in]	stream
      */
      LogStream(LogStreamBuf * buf = nullptr, bool delete_buf = true, std::ostream * stream = nullptr);

      /// Clears all message buffers.
      ~LogStream() override;
      //@}

      /// @name Stream Methods
      //@{

      /**
        rdbuf method of ostream.
        This method is needed to access the LogStreamBuf object.

        @see std::ostream::rdbuf for more details.
      */
      LogStreamBuf * rdbuf();

      /// Arrow operator.
      LogStreamBuf * operator->();
      //@}


      /// @name Level methods
      //@{

      /**
        Set the level of the LogStream

       @param[in] level The new LogLevel
      */
      void setLevel(std::string level);


      /**
        Returns the LogLevel of this LogStream
      */
      std::string getLevel();
      //@}

      /**
        @name Associating Streams

        A thread-local LogStream (see getThreadLocalLog*()) writes to the destinations of the
        corresponding global LogStream, including later changes to them. Its own changes apply on
        top of these, to the calling thread only: a stream inserted into it belongs to it and is not
        affected by later changes to the global LogStream; a global destination removed from it stays
        hidden until it is inserted into it again. A destination that the global LogStream inserts later,
        also the same stream again, is not hidden; a LogSinkGuard on the global LogStream restores its
        destination as it was, so it stays hidden.
        setPrefix() on it applies to its own destinations only; the prefix of a global destination is
        set on the global LogStream.

        Do not insert a LogStream as a destination of another one: writing to it would block.

        When the destinations change, pending repeat counts ("<message> occurred N times") are written to the
        destinations that got the repeated messages and remain (with insert(), remove() and removeAllStreams() on
        this stream: to all destinations before the change), and the earlier messages no longer count as repeated:
        a new destination gets every message from then on.
      */
      //@{

      /**
        Associate a new stream with this logstream.
        This method inserts a new stream into the list of
        associated streams and sets the corresponding minimum
        and maximum log levels.
        Any message that is subsequently logged, will be copied
        to this stream if its log level is between <tt>min_level</tt>
        and <tt>max_level</tt>. If <tt>min_level</tt> and <tt>max_level</tt>
        are omitted, all messages are copied to this stream.
        If <tt>min_level</tt> and <tt>max_level</tt> are equal, this function can be used
        to listen to a specified channel.

        @param[in] s a reference to the stream to be associated
      */
      void insert(std::ostream & s);

      /**
        Remove an association with a stream.

        Remove a stream from the stream list and avoid the copying of new messages to
        this stream. \par
        If the stream was not in the list of associated streams nothing will
        happen.

        @param[in] s the stream to be removed
      */
      void remove(std::ostream & s);

      /**
        Remove all streams associated to this LogStream, effectively silencing it.
        
        Flushes all buffers to ensure any pending log messages are written to their
        respective streams before they are removed.
      */
      void removeAllStreams();

      /**
        Set prefix for output to this stream.
        Each line written to the stream will be prefixed by
        this string. The string may also contain trivial
        format specifiers to include loglevel and time/date
        of the logged message. \par
        The following format tags are recognized:

        - <b>%y</b> message type ("Error", "Warning", "Information", "-")
        - <b>%T</b> time (HH:MM:SS)
        - <b>%t</b> time in short format (HH:MM)
        - <b>%D</b>	date (YYYY/MM/DD)
        - <b>%d</b> date in short format (MM/DD)
        - <b>%S</b> time and date (YYYY/MM/DD, HH:MM:SS)
        - <b>%s</b> time and date in short format (MM/DD, HH:MM)
        - <b>%%</b>	percent sign (escape sequence)

        @param[in] s The stream that will be prefixed.
        @param[in] prefix The prefix used for the stream.
      */
      void setPrefix(const std::ostream & s, const std::string & prefix);


      /// Set prefix of all output streams, details see setPrefix method with ostream
      void setPrefix(const std::string & prefix);

      ///
      void flush();

      /**
        @brief Flush any incomplete line (text not terminated by newline) to all streams.

        This is useful for thread-local logging where the buffer needs to be
        flushed before streams are reconfigured globally.
      */
      void flushIncomplete();

      /// Is @p stream currently one of this LogStream's output destinations?
      bool hasStream(std::ostream & stream);
      //@}
private:

      typedef std::list<LogStreamBuf::StreamStruct>::iterator StreamIterator;

      StreamIterator findStream_(const std::ostream & stream);
      bool hasStream_(std::ostream & stream);
      bool bound_() const;

      /// Lets the thread-local accessors (getThreadLocalLog*()) call followParent_() on each use
      friend struct ThreadLocalLogAccess;

      friend class LogSinkGuard;

      /// A destination that LogSinkGuard removed from a global stream (see removeForGuard_()), restored with its prefix and id
      struct GuardedRemoval
      {
        bool reinsert = false;
        std::string prefix;
        Size id = 0;
      };

      /**
        For LogSinkGuard: on a thread-local stream, suppresses @p s until restoreForGuard_(); on a global stream, removes it.
        Pending text and repeat counts go to the current destinations first; for a global stream also those of the calling
        thread's thread-local stream of the same level, which follows it.
      */
      GuardedRemoval removeForGuard_(std::ostream& s);

      /**
        Undoes removeForGuard_(). Text pending since then goes to the destinations without @p s first. A thread-local stream
        then has the destinations it would have without the guard; a global stream gets @p s back unless it was inserted again.
      */
      void restoreForGuard_(std::ostream& s, const GuardedRemoval& removal);

      /// Removes @p s from the destinations (for a thread-local stream: hides inherited ones) and returns what it removed. The caller holds the sink mutex.
      GuardedRemoval detachLocked_(std::ostream& s);

      /// For a global stream: the calling thread's thread-local stream of the same level, or nullptr
      LogStream* threadLocalFollower_();

      /// Writes pending text (also an incomplete line) and repeat counts to the current destinations
      void drain_();

      /// For a thread-local stream: applies changes of the global destinations made since the last call, see updateState_()
      void followParent_();

      /**
        For a thread-local stream: sets badbit while it has no destination, so that operator<< does not format
        messages that would be discarded anyway, and clears the state once it has a destination again.
        Global streams are left as they are; a thread-local stream is used by its own thread only.
      */
      void updateState_();

      /// flag needed by the destructor to decide whether the streambuf
      /// has to be deleted. If the default ctor is used to create
      /// the LogStreamBuf, delete_buffer_ is set to true and the ctor
      /// also deletes the buffer.
      bool delete_buffer_;

    }; //LogStream

    /**
      @brief RAII guard that temporarily removes an attached stream from a LogStream and re-inserts it on scope exit.

      This class provides exception-safe temporary removal of output streams from LogStream objects.
      Use this when you need to temporarily suppress logging to a specific stream (e.g., cout)
      during operations that may throw exceptions.

      @note The guard acts only on a stream that is attached when it is constructed; for an
            unattached stream it does nothing (will not insert on scope exit). On a thread-local stream,
            it suppresses the stream for the calling thread until it ends, also if the global stream
            inserts it again meanwhile; then the thread-local stream has the destinations it would have
            without the guard (e.g. not a destination that the global stream removed meanwhile). On a
            global stream, it removes the destination and inserts it again, with its prefix, at the end
            (unless it was inserted again meanwhile); the stream must outlive the guard.

      @note The OPENMS_LOG_* macros write to thread-local streams that follow the destinations of
            the global ones. Guarding a global stream (getGlobalLog*()) suppresses the sink for all
            threads; guarding a thread-local stream (getThreadLocalLog*()) only for the calling thread.
            Text pending before and inside the guarded scope goes to the destinations it was written
            for: for a global stream, this holds for the calling thread, not for other threads.

      Example usage:
      @code
      LogSinkGuard guard(getThreadLocalLogInfo(), cout); // Removes cout, will re-insert on scope exit
      some_operation_that_may_throw();
      // cout is automatically re-inserted when guard goes out of scope
      @endcode

      @ingroup Concept
    */
    class OPENMS_DLLAPI LogSinkGuard
    {
    public:
      /**
        @brief Construct a guard that removes @p stream, if attached, and re-inserts it on destruction.

        @param log_stream The LogStream to remove the stream from (and re-insert into on destruction)
        @param stream The stream to temporarily remove (e.g., std::cout); if it is not one of
                      @p log_stream's destinations, the guard does nothing
      */
      LogSinkGuard(LogStream& log_stream, std::ostream& stream)
        : log_stream_(log_stream), stream_(stream), was_attached_(log_stream.hasStream(stream))
      {
        if (!was_attached_)
        {
          return;
        }
        // Drains pending output first, so that pre-guard text reaches this sink before it is detached.
        removal_ = log_stream_.removeForGuard_(stream_);
      }

      /// Destructor writes output pending inside the scope to the other destinations, then restores the sink.
      ~LogSinkGuard()
      {
        if (!was_attached_)
        {
          return;
        }
        log_stream_.restoreForGuard_(stream_, removal_);
      }

      // Non-copyable and non-movable
      LogSinkGuard(const LogSinkGuard&) = delete;
      LogSinkGuard& operator=(const LogSinkGuard&) = delete;
      LogSinkGuard(LogSinkGuard&&) = delete;
      LogSinkGuard& operator=(LogSinkGuard&&) = delete;

    private:
      LogStream& log_stream_;
      std::ostream& stream_;
      /// was the sink attached when the guard was constructed? if not, the guard does nothing
      const bool was_attached_;
      /// what the guard removed, restored on destruction
      LogStream::GuardedRemoval removal_;
    };

  } // namespace Logger

  //
  // Thread-Local Log Stream Accessors
  //
  // These functions return thread-local LogStream instances, eliminating data races
  // in the logging system. Each thread gets its own buffers and state.
  // See GitHub Issue #8596 for details on the race conditions this fixes.
  //
  // Each thread-local stream writes to the destinations of the corresponding global stream
  // (getGlobalLog*()). A change of the global destinations (e.g. via LogConfigHandler) takes
  // effect in every thread at its next message, including threads that have logged before;
  // a destination removed from the global stream is not written to afterwards, so it may be
  // destroyed. Changes of a thread-local stream itself apply on top, to its thread only
  // (see LogStream::insert()).
  //
  // A thread-local stream without destinations (e.g. the debug stream, unless debug output is enabled) does
  // not format messages: it is in a failed state (badbit) until it has a destination again. The arguments of
  // a message are still evaluated, except with OPENMS_LOG_DEBUG, which skips the whole message then.
  // Each call of these accessors (as in the OPENMS_LOG_* macros) applies changes
  // of the global destinations to the calling thread's stream; a reference kept across such a change picks it up
  // with the next call on that thread, or with insert() or remove() on it.
  //
  // RESTRICTIONS:
  // - Change the global configuration from one thread at a time.
  // - OpenMP thread pools reuse threads, so thread_local state (e.g. local changes) persists across parallel regions.
  //

  /// @brief Get thread-local fatal error log stream
  OPENMS_DLLAPI Logger::LogStream& getThreadLocalLogFatal();

  /// @brief Get thread-local error log stream
  OPENMS_DLLAPI Logger::LogStream& getThreadLocalLogError();

  /// @brief Get thread-local warning log stream
  OPENMS_DLLAPI Logger::LogStream& getThreadLocalLogWarn();

  /// @brief Get thread-local info log stream
  OPENMS_DLLAPI Logger::LogStream& getThreadLocalLogInfo();

  /// @brief Get thread-local debug log stream
  OPENMS_DLLAPI Logger::LogStream& getThreadLocalLogDebug();

  /**
    @brief Enables or disables writing OPENMS_LOG_DEBUG messages to std::cout.

    Attaches std::cout to, or detaches it from, the global debug stream, which all threads follow.
    It also undoes a local removal (or insertion) of std::cout on the calling thread's debug stream.
    Pending and cached messages of the calling thread are flushed to the previous destinations first.

    @param enabled Attach std::cout if true, detach it otherwise
  */
  OPENMS_DLLAPI void setConsoleDebugLogging(bool enabled);

  //
  // Logging Macros (thread-safe via thread-local streams)
  //

  /// Macro for fatal errors (processing stops) - includes file and line info
#define OPENMS_LOG_FATAL_ERROR \
  OpenMS::getThreadLocalLogFatal() << __FILE__ << "(" << __LINE__ << "): "

  /// Macro for non-fatal errors (processing continues)
#define OPENMS_LOG_ERROR \
  OpenMS::getThreadLocalLogError()

  /// Macro for warnings
#define OPENMS_LOG_WARN \
  OpenMS::getThreadLocalLogWarn()

  /// Macro for information/status messages
#define OPENMS_LOG_INFO \
  OpenMS::getThreadLocalLogInfo()

  namespace Logger
  {
    /// Turns a whole message into a void expression, which OPENMS_LOG_DEBUG skips when debug output is disabled
    struct LogVoidify
    {
      void operator&(std::ostream&) const {}
    };
  }

  /**
    Macro for debug information - includes file and line info.

    Without a destination for debug output on the calling thread (the default, unless e.g. a TOPP tool runs
    with -debug), the whole message is skipped: its arguments are not evaluated. As it is a void expression,
    use getThreadLocalLogDebug() for anything else than writing a message.
  */
#define OPENMS_LOG_DEBUG \
  !OpenMS::getThreadLocalLogDebug().good() ? (void)0 : OpenMS::Logger::LogVoidify() & \
  OpenMS::getThreadLocalLogDebug() << past_last_slash(__FILE__) << "(" << __LINE__ << "): "

  /// Macro for debug information (without file info), see OPENMS_LOG_DEBUG
#define OPENMS_LOG_DEBUG_NOFILE \
  !OpenMS::getThreadLocalLogDebug().good() ? (void)0 : OpenMS::Logger::LogVoidify() & OpenMS::getThreadLocalLogDebug()

  /**
    @name Global LogStream accessor functions

    These functions provide access to global LogStream instances for configuration purposes
    (e.g., adding/removing output streams). For actual logging, use the OPENMS_LOG_* macros
    which use thread-local streams to avoid data races.

    @warning Direct logging to global streams is NOT thread-safe. Use OPENMS_LOG_* macros instead.
  */
  ///@{

  /**
    @brief Get the global fatal error log stream for configuration purposes.

    Returns a reference to the global fatal error LogStream instance. Use this function
    to configure the stream (e.g., adding/removing output destinations with insert()/remove()).

    @warning Do NOT use this for logging messages directly - it is not thread-safe.
             Use the OPENMS_LOG_FATAL_ERROR macro instead.

    @return Reference to the global fatal error log stream

    @see OPENMS_LOG_FATAL_ERROR
  */
  OPENMS_DLLAPI Logger::LogStream& getGlobalLogFatal();

  /**
    @brief Get the global error log stream for configuration purposes.

    Returns a reference to the global error LogStream instance. Use this function
    to configure the stream (e.g., adding/removing output destinations with insert()/remove()).

    @warning Do NOT use this for logging messages directly - it is not thread-safe.
             Use the OPENMS_LOG_ERROR macro instead.

    @return Reference to the global error log stream

    @see OPENMS_LOG_ERROR
  */
  OPENMS_DLLAPI Logger::LogStream& getGlobalLogError();

  /**
    @brief Get the global warning log stream for configuration purposes.

    Returns a reference to the global warning LogStream instance. Use this function
    to configure the stream (e.g., adding/removing output destinations with insert()/remove()).

    @warning Do NOT use this for logging messages directly - it is not thread-safe.
             Use the OPENMS_LOG_WARN macro instead.

    @return Reference to the global warning log stream

    @see OPENMS_LOG_WARN
  */
  OPENMS_DLLAPI Logger::LogStream& getGlobalLogWarn();

  /**
    @brief Get the global info log stream for configuration purposes.

    Returns a reference to the global info LogStream instance. Use this function
    to configure the stream (e.g., adding/removing output destinations with insert()/remove()).

    @warning Do NOT use this for logging messages directly - it is not thread-safe.
             Use the OPENMS_LOG_INFO macro instead.

    @return Reference to the global info log stream

    @see OPENMS_LOG_INFO
  */
  OPENMS_DLLAPI Logger::LogStream& getGlobalLogInfo();

  /**
    @brief Get the global debug log stream for configuration purposes.

    Returns a reference to the global debug LogStream instance. Use this function
    to configure the stream (e.g., adding/removing output destinations with insert()/remove()).

    @note The debug stream is disabled by default. It is enabled in TOPPBase when running
          with '-debug' (see setConsoleDebugLogging()).

    @warning Do NOT use this for logging messages directly - it is not thread-safe.
             Use the OPENMS_LOG_DEBUG macro instead.

    @return Reference to the global debug log stream

    @see OPENMS_LOG_DEBUG
  */
  OPENMS_DLLAPI Logger::LogStream& getGlobalLogDebug();

  ///@}

} // namespace OpenMS
