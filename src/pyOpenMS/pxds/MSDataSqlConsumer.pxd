from Types cimport *
from MSChromatogram cimport *
from MSSpectrum cimport *
from ExperimentalSettings cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>" namespace "OpenMS":
    
    cdef cppclass MSDataSqlConsumer:

        MSDataSqlConsumer(String filename, UInt64 run_id, int buffer_size, bool full_meta, bool lossy_compression, double linear_mass_acc) except + nogil 
        MSDataSqlConsumer(MSDataSqlConsumer &) except + nogil  # compiler

        void flush() except + nogil 
            # wrap-doc:
                #  Flushes the data for good
                #  
                #  After calling this function, no more data is held in the buffer but the
                #  class is still able to receive new data

        void consumeSpectrum(MSSpectrum & s) except + nogil  # wrap-doc:Write a spectrum to the output file
        void consumeChromatogram(MSChromatogram & c) except + nogil  # wrap-doc:Write a chromatogram to the output file

        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) except + nogil 
        void setExperimentalSettings(ExperimentalSettings & exp) except + nogil 
