from SourceFile cimport *
from DateTime cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/METADATA/ExperimentalSettings.h>" namespace "OpenMS":

    cdef cppclass ExperimentalSettings:

        ExperimentalSettings() nogil except +
        ExperimentalSettings(ExperimentalSettings) nogil except + # wrap-ignore

        libcpp_vector[SourceFile] getSourceFiles() nogil except +
        void setSourceFiles(libcpp_vector[SourceFile] source_files) nogil except +

        DateTime getDateTime() nogil except +
        void setDateTime(DateTime date_time) nogil except +

        String getLoadedFilePath() nogil except +
        void   setLoadedFilePath(String from_) nogil except +
