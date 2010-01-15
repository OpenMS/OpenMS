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

#ifndef OPENMS_CONCEPT_LOGSTREAM_H
#define OPENMS_CONCEPT_LOGSTREAM_H

#include <OpenMS/DATASTRUCTURES/String.h>


#include <sstream>
#include <iostream>
#include <list>
#include <vector>
#include <ctime>
#include <map>

namespace OpenMS 
{
	namespace Logger
	{

	/**	@name Log streams
			Logging, filtering, and storing messages.
			Many programs emit warning messages, error messages, or simply
			informations and remarks to their users. The  LogStream  
			class provides a convenient and straight-forward interface 
			to classify these messages according to their importance 
			(via the loglevel), filter and store them in files or
			write them to streams. \par
			As the LogStream class is derived from ostream, it behaves 
			as any ostream object. Additionally you may associate
			streams with each LogStream object that catch only 
			messages of certain loglevels. So the user might decide to
			redirect all error messages to cerr, all warning messages
			to cout and all informations to a file. \par
			Along with each message its time of creation and its loglevel
			is stored. So the user might also decide to store all 
			errors he got in the last two hours or alike. \par
			The LogStream class heavily relies on the LogStreamBuf 
			class, which does the actual buffering and storing, but is only
			of interest if you want to implement a derived class, as the 
			actual user interface is implemented in the LogStream class.
	* 	\ingroup Common
	*/
	//@{

	// forward declarations
	class LogStream;
	class LogStreamNotifier;

	/** Stream buffer used by LogStream.
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
	class OPENMS_DLLAPI LogStreamBuf
		: public std::streambuf
	{

		friend class LogStream;

		public:

		/**	@name	Constants
		*/
		//@{
		static const time_t MAX_TIME;
		static const std::string UNKNOWN_LOG_LEVEL;
		//@}

		/**	@name Constructors and Destructors
		*/
		//@{
		
		/** Default constructor.
				Create a new LogStreamBuf object. The level is set to unknown
		*/
		LogStreamBuf();

    /**
        Create a new LogStreamBuf object and set the level to log_level
        @param log_level The log level of the LogStreamBuf
    */
    LogStreamBuf(std::string log_level);

		/** Destructor.
				Destruct the buffer and free all stored messages strings.
		*/
		virtual ~LogStreamBuf();
		
		//@}
		
		/**	@name	Stream methods 
		*/
		//@{

		/**	Sync method.
				This method is called as soon as the ostream is flushed
				(especially this method is called by flush or endl).
				It transfers the contents of the streambufs putbuffer 
				into a logline if a newline or linefeed character 
				is found in the buffer ("\n" or "\r" resp.).
				The line is then removed from the putbuffer.
				Incomplete lines (not terminated by "\n" / "\r" are
				stored in incomplete_line_.
		*/
		virtual int sync();

		/**	Overflow method.
				This method calls sync and <tt>streambuf::overflow(c)</tt> to 
				prevent a buffer overflow.
		*/
		virtual int overflow(int c = -1);
		//@}


    /**
     * @name Level methods
     */
    //@{

    /**
     * Set the level of the LogStream
     *
     * @param level The new LogLevel
     */
    void setLevel(std::string level);


    /**
     * Returns the LogLevel of this LogStream
     */
    std::string getLevel();
    //@}

    /**
     * @brief Holds a stream that is connected to the LogStream incl. the minimum and maximum
     * level at which the LogStream redirects messages to this stream.
     */
		OPENMS_DLLAPI struct StreamStruct
		{
			std::ostream*				stream;
      std::string         prefix;
			LogStreamNotifier*	target;
		
			StreamStruct()
				:	stream(0),
					target(0)
			{
			}
			
			// Delete the notification target.
			~StreamStruct()
			{
			}
		};


		protected:
		
		/// distribute a new message to connected streams
		void distribute_(std::string outstring);
		
		// interpret the prefix format string and return the expanded prefix
		std::string expandPrefix_(const std::string& prefix, time_t time) const;

		char* 									pbuf_;
		std::string             level_;
		std::list<StreamStruct>	stream_list_;
		std::string             incomplete_line_;

 		/**	@name Caching
		*/
		//@{

    /**
     * @brief Holds a counter of occurences and an index for the occurence sequence 
     * of the corresponding log message
     */
    struct LogCacheStruct 
    {
      Size timestamp;
      int counter;
    };

    /**
     * Sequential counter to remember the sequence of occurence 
     * of the cached log messages
     */
    Size log_cache_counter_;

    /// Cache of the last two log messages
    std::map<std::string, LogCacheStruct> log_cache_;
    /// Cache of the occurence sequence of the last two log messages
    std::map<Size, std::string > log_time_cache_;

    /// Checks if the line is already in the cache
    bool isInCache_(std::string const & line);

    /**
       Adds the new line to the cache and removes an old one
       if necessary

       @param line The Log message that should be added to the cache
       @return An additional massage if a reoccuring message was removed
       from the cache
     */
    std::string addToCache_(std::string const & line);

    /**
     * Returns the next free index for a log message
     */
    Size getNextLogCounter_();

    /**
     * Checks if some of the cached entries where sent more then once
     * to the LogStream and (if necessary) prints a corresponding messages
     * into all affected Logs
     */
    void clearCache_();
		//@}

	};


	///
	class OPENMS_DLLAPI LogStreamNotifier
	{
		public:
		
		///
		LogStreamNotifier();
			
		///
		virtual ~LogStreamNotifier();

		///
		virtual void logNotify();

		///
		void registerAt(LogStream& log_stream);
		///
		void unregister();

		protected:

		std::stringstream stream_;

		LogStream* registered_at_;
	};



	/**	Log Stream Class.

	Defines a log stream which features a cache and some formatting.
	For the developer, however, only some macros are of interest which
	will push the message that follows them into the
	appropriate stream:
	
	Macros:
		- LOG_FATAL_ERROR
		- LOG_ERROR (non-fatal error are reported (processing continues))
		- LOG_WARN  (warning, a piece of information which should be read by the user, should be logged)
		- LOG_INFO (information, e.g. a status should be reported)
		- LOG_DEBUG (general debugging information)
				 
			To use a specific logger of a log level simply
			use it as cerr or cout: <br>
			<code> LOG_ERROR << " A bad error occured ..."  </code> <br>
			Which produces an error message in the log.
	*/
			 

	class OPENMS_DLLAPI LogStream
		: public std::ostream
	{
		public:

		/**	@name	Constructors and Destructors
		*/
		//@{

		/** Constructor.
				Creates a new LogStream object that is not associated with any stream.
				If the argument <tt>stream</tt> is set to an output stream (e.g. <tt>cout</tt>)
				all output is send to that stream.
				@param	buf
				@param  delete_buf
				@param	stream
		*/
		LogStream(LogStreamBuf* buf = 0, bool delete_buf = true, std::ostream* stream = 0);

		/** Destructor.
				Clears all message buffers.
		*/
		virtual ~LogStream();
	
		//@}		

		/**	@name	Stream Methods
		*/
		//@{

		/**	<tt>rdbuf</tt> method of ostream.
				This method is needed to access the LogStreamBuf object.
		*/
		LogStreamBuf* rdbuf();

		/** Arrow operator.
		*/
		LogStreamBuf* operator -> ();
		//@}


		/**
		 * @name Level methods
		 */
		//@{

		/**
		 * Set the level of the LogStream
		 *
		 * @param level The new LogLevel
		 */
		void setLevel(std::string level);


		/**
		 * Returns the LogLevel of this LogStream
		 */
		std::string getLevel();
		//@}

		/**	@name Associating Streams 
		*/
		//@{

		/**	Associate a new stream with this logstream.
				This method inserts a new stream into the list of 
				associated streams and sets the corresponding minimum
				and maximum log levels.
				Any message that is subsequently logged, will be copied
				to this stream if its log level is between <tt>min_level</tt>
				and <tt>max_level</tt>. If <tt>min_level</tt> and <tt>max_level</tt>
				are omitted, all messages are copied to this stream.
				If <tt>min_level</tt>	and <tt>max_level</tt> are equal, this function can be used
				to listen to a specified channel.
				@param	s a reference to the stream to be associated
		*/
		void insert
			(std::ostream& s);

		/**	Remove an association with a stream.
				Remove a stream from the stream list and avoid the copying of new messages to
				this stream. \par
				If the stream was not in the list of associated streams nothing will
				happen.
				@param	s the stream to be removed
		*/
		void remove(std::ostream& s);

		/**	Add a notification target
		*/
		void insertNotification(std::ostream& s, 
														LogStreamNotifier& target);

		/**	Set prefix for output to this stream.
				Each line written to the stream will be prefixed by
				this string. The string may also contain trivial 
				format specifiers to include loglevel and time/date 
				of the logged message. \par
				The following format tags are recognized:

					- <b>%y</b>	message type ("Error", "Warning", "Information", "-")
					- <b>%T</b> time (HH:MM:SS)
					- <b>%t</b>  time in short format (HH:MM)
					- <b>%D</b>	date (DD.MM.YYYY)
					- <b>%d</b>  date in short format (DD.MM.)
					- <b>%S</b> time and date (DD.MM.YYYY, HH:MM:SS)
					- <b>%s</b>  time and date in short format (DD.MM., HH:MM)
					- <b>%%</b>	percent sign (escape sequence)
				
		*/
		void setPrefix(const std::ostream& s, const std::string& prefix);
		

		///	Set prefix of all output streams, details see setPrefix method with ostream
		void setPrefix(const std::string& prefix);

		///
		void flush();
		//@}		

		private:

		typedef std::list<LogStreamBuf::StreamStruct>::iterator StreamIterator;
		
		StreamIterator findStream_(const std::ostream& stream);
		bool hasStream_(std::ostream& stream);
		bool bound_() const;

		/// flag needed by the destructor to decide whether the streambuf
		/// has to be deleted. If the default ctor is used to create
		/// the LogStreamBuf, delete_buffer_ is set to true and the ctor
		/// also deletes the buffer.
		bool	delete_buffer_;

	}; //LogStream

	} // namespace Logger

	
	/// Macro to be used if fatal error are reported (processing stops)
	#define LOG_FATAL_ERROR \
  Log_fatal << __FILE__ << "(" << __LINE__ << "): "
	
	/// Macro to be used if non-fatal error are reported (processing continues)
	#define LOG_ERROR \
  Log_error

	/// Macro if a warning, a piece of information which should be read by the user, should be logged
  #define LOG_WARN \
  Log_warn

	/// Macro if a information, e.g. a status should be reported
  #define LOG_INFO \
  Log_info

	/// Macro for general debugging information
  #define LOG_DEBUG \
  Log_debug << __FILE__ << "(" << __LINE__ << "): "


	/**
 Global static instance of a logstream.
			This instance of LogStream is by default bound to <b>cout</b> <b>cerr</b> by calling
			the default constructor.
	*/

	OPENMS_DLLAPI extern Logger::LogStream	Log_fatal;
	OPENMS_DLLAPI extern Logger::LogStream  Log_error;
	OPENMS_DLLAPI extern Logger::LogStream  Log_warn;
	OPENMS_DLLAPI extern Logger::LogStream  Log_info;
	OPENMS_DLLAPI extern Logger::LogStream  Log_debug;

} // namespace OpenMS

#endif // OPENMS_CONCEPT_LOGSTREAM_H
