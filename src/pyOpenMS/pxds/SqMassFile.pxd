from Types cimport *
from libcpp cimport bool
from MzMLSqliteHandler cimport *
from Types cimport *
from MSExperiment cimport *
from IMSDataConsumer cimport *

cdef extern from "<OpenMS/FORMAT/SqMassFile.h>" namespace "OpenMS":
    
    cdef cppclass SqMassFile:
        SqMassFile() nogil except +
        SqMassFile(SqMassFile) nogil except + #wrap-ignore
        void load(const String & filename, MSExperiment & map_) nogil except +
        void store(const String & filename, MSExperiment & map_) nogil except +
        # NAMESPACE # # POINTER # void transform(const String & filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count, bool skip_first_pass) nogil except +
        void setConfig(SqMassConfig config) nogil except +

cdef extern from "<OpenMS/FORMAT/SqMassFile.h>" namespace "OpenMS::SqMassFile":
    
    cdef cppclass SqMassConfig "OpenMS::SqMassFile::SqMassConfig":
        SqMassConfig() nogil except +
        SqMassConfig(SqMassConfig) nogil except + #wrap-ignore
        bool write_full_meta
        bool use_lossy_numpress
        double linear_fp_mass_acc

