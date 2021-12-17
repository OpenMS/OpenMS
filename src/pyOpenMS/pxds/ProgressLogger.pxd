from String cimport *
from Types cimport *

cdef extern from "<OpenMS/CONCEPT/ProgressLogger.h>" namespace "OpenMS":

    cdef cppclass ProgressLogger:
        # wrap-doc:
            #   Base class for all classes that want to report their progress
            #   -----
            #   Per default the progress log is disabled. Use setLogType to enable it
            #   -----
            #   Use startProgress, setProgress and endProgress for the actual logging

        ProgressLogger() nogil except + 
        ProgressLogger(ProgressLogger &) nogil except +
        void setLogType(LogType) nogil except + # wrap-doc:Sets the progress log that should be used. The default type is NONE!
        LogType getLogType() nogil except + # wrap-doc:Returns the type of progress log being used
        void startProgress(SignedSize begin, SignedSize end, String label) nogil except +
        void setProgress(SignedSize value) nogil except + # wrap-doc:Sets the current progress
        void endProgress() nogil except + # wrap-doc:Ends the progress display
        void nextProgress() nogil except + # wrap-doc:Increment progress by 1 (according to range begin-end)

cdef extern from "<OpenMS/CONCEPT/ProgressLogger.h>" namespace "OpenMS::ProgressLogger":

    cdef enum LogType:
        # wrap-attach: ProgressLogger
        CMD, GUI, NONE
