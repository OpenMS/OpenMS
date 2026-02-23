from libcpp.vector cimport vector

from QCBase cimport QCBase
from MSExperiment cimport MSExperiment
from MSChromatogram cimport MSChromatogram
from MzTabMetaData cimport MzTabMetaData
from String cimport String
from Types cimport UInt


cdef extern from "<OpenMS/QC/TIC.h>" namespace "OpenMS":

    cdef cppclass TIC(QCBase):
        TIC() except +
        TIC_Result compute(MSExperiment& exp, float bin_size, unsigned int ms_level) except +
        const String& getName() const

    cdef cppclass TIC_Result "OpenMS::TIC::Result":
        vector[unsigned int] intensities
        vector[float] relative_intensities
        vector[float] retention_times
        unsigned int area
        unsigned int fall
        unsigned int jump

        TIC() except +

        Result compute(const MSExperiment& exp,
                       float bin_size = 0,
                       UInt ms_level = 1) except +

        const String& getName() const except +

        const vector[MSChromatogram]& getResults() const except +

        QCBase.Status requirements() const except +

        void addMetaDataMetricsToMzTab(MzTabMetaData& meta,
                                       vector[Result]& tics) except +
