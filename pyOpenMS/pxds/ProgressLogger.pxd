from String cimport *
from Types cimport *

cdef extern from "<OpenMS/CONCEPT/ProgressLogger.h>" namespace "OpenMS":

    cdef cppclass ProgressLogger:
        ProgressLogger()           nogil except +
        void setLogType(LogType)           nogil except +
        LogType getLogType()           nogil except +
        void startProgress(SignedSize begin, SignedSize end, String label)           nogil except +
        void setProgress(SignedSize value)           nogil except +
        void endProgress()           nogil except +


cdef extern from "<OpenMS/CONCEPT/ProgressLogger.h>" namespace "OpenMS::ProgressLogger":

    cdef enum LogType:
        # wrap-attach: ProgressLogger
        CMD, GUI, NONE


