from libcpp.vector cimport vector

from QCBase cimport QCBase
from MSExperiment cimport MSExperiment
from MSChromatogram cimport MSChromatogram
from MzTabMetaData cimport MzTabMetaData
from String cimport String
from Types cimport UInt


cdef extern from "<OpenMS/QC/TIC.h>" namespace "OpenMS":

    #  Nested struct first 
    cdef cppclass TIC_Result "OpenMS::TIC::Result":
        vector[UInt] intensities
        vector[float] relative_intensities
        vector[float] retention_times
        UInt area
        UInt fall
        UInt jump


    
    cdef cppclass TIC(QCBase):

        TIC() except +

        TIC_Result compute(const MSExperiment& exp,
                           float bin_size = 0,
                           UInt ms_level = 1) except +

        const String& getName() const

        const vector[MSChromatogram]& getResults() except + const

        QCBase.Status requirements() except + const

        void addMetaDataMetricsToMzTab(MzTabMetaData& meta,
                                       vector[TIC_Result]& tics) except +
