// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Stephan Aiche, Andreas Bertsch$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/String.h>

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
    class LogStreamNotifier;

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

        @param log_level The log level of the LogStreamBuf (default is unknown)
        @param col If messages should be colored, provide a colorizer here
      */
      LogStreamBuf(const std::string& log_level = UNKNOWN_LOG_LEVEL, Colorizer* col = nullptr);

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

        @param level The new LogLevel
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
        LogStreamNotifier * target;

        StreamStruct() :
          stream(nullptr),
          target(nullptr)
        {}

        /// Delete the notification target.
        ~StreamStruct()
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

      /// Interpret the prefix format string and return the expanded prefix.
      std::string expandPrefix_(const std::string & prefix, time_t time) const;

      char * pbuf_ = nullptr;
      std::string             level_;
      std::list<StreamStruct> stream_list_;
      std::string             incomplete_line_;
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

        @param line The Log message that should be added to the cache
        @return An additional massage if a re-occurring message was removed
        from the cache
      */
      std::string addToCache_(std::string const & line);

      /// Returns the next free index for a log message
      Size getNextLogCounter_();

      /// Non-lock acquiring sync function called in the d'tor
      int syncLF_();
      //@}
    };

    ///
    class OPENMS_DLLAPI LogStreamNotifier
    {
public:

      /// Empty constructor.
      LogStreamNotifier();

      /// Destructor
      virtual ~LogStreamNotifier();

      ///
      virtual void logNotify();

      ///
      void registerAt(LogStream & log_stream);

      ///
      void unregister();

protected:
      std::stringstream stream_;

      LogStream * registered_at_;
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

      @note The log stream macros are thread safe and can be used in a
      multithreaded environment, the global variables are not! The macros are
      protected by a OPENMS_THREAD_CRITICAL directive (which translates to an
      OpenMP critical pragma), however there may be a small performance penalty
      to this.

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

        @param	buf
        @param  delete_buf
        @param	stream
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

       @param level The new LogLevel
      */
      void setLevel(std::string level);


      /**
        Returns the LogLevel of this LogStream
      */
      std::string getLevel();
      //@}

      /// @name Associating Streams
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

        @param s a reference to the stream to be associated
      */
      void insert(std::ostream & s);

      /**
        Remove an association with a stream.

        Remove a stream from the stream list and avoid the copying of new messages to
        this stream. \par
        If the stream was not in the list of associated streams nothing will
        happen.

        @param s the stream to be removed
      */
      void remove(std::ostream & s);

      /**
        Remove all streams associated to this LogStream, effectively silencing it.
      */
      void removeAllStreams();

      /// Add a notification target
      void insertNotification(std::ostream & s,
                              LogStreamNotifier & target);

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

        @param s The stream that will be prefixed.
        @param prefix The prefix used for the stream.
      */
      void setPrefix(const std::ostream & s, const std::string & prefix);


      /// Set prefix of all output streams, details see setPrefix method with ostream
      void setPrefix(const std::string & prefix);

      ///
      void flush();
      //@}
private:

      typedef std::list<LogStreamBuf::StreamStruct>::iterator StreamIterator;

      StreamIterator findStream_(const std::ostream & stream);
      bool hasStream_(std::ostream & stream);
      bool bound_() const;

      /// flag needed by the destructor to decide whether the streambuf
      /// has to be deleted. If the default ctor is used to create
      /// the LogStreamBuf, delete_buffer_ is set to true and the ctor
      /// also deletes the buffer.
      bool delete_buffer_;

    }; //LogStream

  } // namespace Logger

  /// Macro to be used if fatal error are reported (processing stops)
#define OPENMS_LOG_FATAL_ERROR \
  OPENMS_THREAD_CRITICAL(LOGSTREAM) \
  OpenMS_Log_fatal << __FILE__ << "(" << __LINE__ << "): "

  /// Macro to be used if non-fatal error are reported (processing continues)
#define OPENMS_LOG_ERROR \
  OPENMS_THREAD_CRITICAL(LOGSTREAM) \
  OpenMS_Log_error

  /// Macro if a warning, a piece of information which should be read by the user, should be logged
#define OPENMS_LOG_WARN \
  OPENMS_THREAD_CRITICAL(LOGSTREAM) \
  OpenMS_Log_warn

  /// Macro if a information, e.g. a status should be reported
#define OPENMS_LOG_INFO \
  OPENMS_THREAD_CRITICAL(LOGSTREAM) \
  OpenMS_Log_info

  /// Macro for general debugging information
#define OPENMS_LOG_DEBUG \
  OPENMS_THREAD_CRITICAL(LOGSTREAM) \
  OpenMS_Log_debug << past_last_slash(__FILE__) << "(" << __LINE__ << "): "

  /// Macro for general debugging information (without information on file)
#define OPENMS_LOG_DEBUG_NOFILE \
  OPENMS_THREAD_CRITICAL(LOGSTREAM) \
  OpenMS_Log_debug

  OPENMS_DLLAPI extern Logger::LogStream OpenMS_Log_fatal; ///< Global static instance of a LogStream to capture messages classified as fatal errors. By default it is bound to @b cerr.
  OPENMS_DLLAPI extern Logger::LogStream OpenMS_Log_error; ///< Global static instance of a LogStream to capture messages classified as errors. By default it is bound to @b cerr.
  OPENMS_DLLAPI extern Logger::LogStream OpenMS_Log_warn;  ///< Global static instance of a LogStream to capture messages classified as warnings. By default it is bound to @b cout.
  OPENMS_DLLAPI extern Logger::LogStream OpenMS_Log_info;  ///< Global static instance of a LogStream to capture messages classified as information. By default it is bound to @b cout.
  OPENMS_DLLAPI extern Logger::LogStream OpenMS_Log_debug; ///< Global static instance of a LogStream to capture messages classified as debug output. By default it is not bound to any output stream. TOPP(AS)Base will connect cout, iff 0 < debug-level

} // namespace OpenMS
