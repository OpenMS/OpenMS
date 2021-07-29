from Types cimport *
from MSChromatogram cimport *
from MSSpectrum cimport *
from ExperimentalSettings cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>" namespace "OpenMS":
    
    cdef cppclass MSDataSqlConsumer:

        MSDataSqlConsumer(MSDataSqlConsumer) nogil except + #wrap-ignore
        MSDataSqlConsumer(String filename, UInt64 run_id, int buffer_size, bool full_meta, bool lossy_compression, double linear_mass_acc) nogil except +

        void flush() nogil except +
            # wrap-doc:
                #   Flushes the data for good
                #   -----
                #   After calling this function, no more data is held in the buffer but the
                #   class is still able to receive new data

        void consumeSpectrum(MSSpectrum & s) nogil except + # wrap-doc:Write a spectrum to the output file
        void consumeChromatogram(MSChromatogram & c) nogil except + # wrap-doc:Write a chromatogram to the output file

        void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) nogil except +
        void setExperimentalSettings(ExperimentalSettings & exp) nogil except +

