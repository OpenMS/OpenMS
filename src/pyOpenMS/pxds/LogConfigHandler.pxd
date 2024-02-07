from Param cimport *
from String cimport *
from StringList cimport *
from Types cimport *

from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/CONCEPT/LogConfigHandler.h>" namespace "OpenMS":
    
    cdef cppclass LogConfigHandler "OpenMS::LogConfigHandler":
    # wrap-manual-memory:
    #  cdef AutowrapPtrHolder[_LogConfigHandler] inst

        # private
        LogConfigHandler() except + nogil  # wrap-ignore
        # private
        LogConfigHandler(LogConfigHandler) except + nogil  # wrap-ignore

        Param parse(const StringList & setting) except + nogil 
        #wrap-doc:
        #  Translates the given list of parameter settings into a LogStream configuration
        #  
        #  Translates the given list of parameter settings into a LogStream configuration.
        #  Usually this list stems from a command line call.
        #  
        #  Each element in the stringlist should follow this naming convention
        #  
        #  <LOG_NAME> <ACTION> <PARAMETER>
        #  
        #  with
        #  - LOG_NAME: DEBUG,INFO,WARNING,ERROR,FATAL_ERROR
        #  - ACTION: add,remove,clear
        #  - PARAMETER: for 'add'/'remove' it is the stream name (cout, cerr or a filename), 'clear' does not require any further parameter
        #  
        #  Example:
        #  `DEBUG add debug.log`
        #  
        #  This function will **not** apply to settings to the log handlers. Use configure() for that.
        #  
        #  :param setting: StringList containing the configuration options
        #  :raises ParseError: In case of an invalid configuration.
        #  :return: Param object containing all settings, that can be applied using the LogConfigHandler.configure() method

        void configure(const Param & param) except + nogil 
        # wrap-doc:
        #  Applies the given parameters (@p param) to the current configuration
        #  
        #  <LOG_NAME> <ACTION> <PARAMETER> <STREAMTYPE>
        #  
        #  LOG_NAME: DEBUG, INFO, WARNING, ERROR, FATAL_ERROR
        #  ACTION: add, remove, clear
        #  PARAMETER: for 'add'/'remove' it is the stream name ('cout', 'cerr' or a filename), 'clear' does not require any further parameter
        #  STREAMTYPE: FILE, STRING (for a StringStream, which you can grab by this name using getStream() )
        #  
        #  You cannot specify a file named "cout" or "cerr" even if you specify streamtype 'FILE' - the handler will mistake this for the
        #  internal streams, but you can use "./cout" to print to a file named cout.
        #  
        #  A classical configuration would contain a list of settings e.g.
        #  
        #  `DEBUG add debug.log FILE`
        #  `INFO remove cout FILE` (FILE will be ignored)
        #  `INFO add string_stream1 STRING`
        #  
        #  :raises ElementNotFound: If the LogStream (first argument) does not exist.
        #  :raises FileNotWritable: If a file (or stream) should be opened as log file (or stream) that is not accessible.
        #  :raises IllegalArgument: If a stream should be registered, that was already registered with a different type.


        void setLogLevel(const String & log_level) except + nogil 
        # wrap-doc:
        #  Sets a minimum log_level by removing all streams from loggers lower than that level.
        #  Valid levels are from low to high: "DEBUG", "INFO", "WARNING", "ERROR", "FATAL_ERROR"


## wrap static methods
cdef extern from "<OpenMS/CONCEPT/LogConfigHandler.h>" namespace "OpenMS::LogConfigHandler":
    
    LogConfigHandler* getInstance() except + nogil  # wrap-ignore
