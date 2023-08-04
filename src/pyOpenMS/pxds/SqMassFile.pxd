from Types cimport *
from libcpp cimport bool
from MzMLSqliteHandler cimport *
from Types cimport *
from MSExperiment cimport *
from IMSDataConsumer cimport *

cdef extern from "<OpenMS/FORMAT/SqMassFile.h>" namespace "OpenMS":
    
    cdef cppclass SqMassFile:
    # wrap-doc:
            #  An class that uses on-disk SQLite database to read and write spectra and chromatograms
            #  
            #  This class provides functions to read and write spectra and chromatograms
            #  to disk using a SQLite database and store them in sqMass format. This
            #  allows users to access, select and filter spectra and chromatograms
            #  on-demand even in a large collection of data

        SqMassFile() except + nogil 
        SqMassFile(SqMassFile &) except + nogil  # compiler
        void load(const String & filename, MSExperiment & map_) except + nogil  # wrap-doc:Read / Write a complete mass spectrometric experiment
        void store(const String & filename, MSExperiment & map_) except + nogil  # wrap-doc:Store an MSExperiment in sqMass format
        # NAMESPACE # # POINTER # void transform(const String & filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count, bool skip_first_pass) except + nogil 
        void setConfig(SqMassConfig config) except + nogil 

cdef extern from "<OpenMS/FORMAT/SqMassFile.h>" namespace "OpenMS::SqMassFile":
    
    cdef cppclass SqMassConfig "OpenMS::SqMassFile::SqMassConfig":
        SqMassConfig() except + nogil  # compiler
        SqMassConfig(SqMassConfig &) except + nogil  # compiler
        bool write_full_meta
        bool use_lossy_numpress
        double linear_fp_mass_acc
