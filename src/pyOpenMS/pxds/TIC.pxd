from libcpp.vector cimport vector

from QCBase cimport QCBase
from MSExperiment cimport MSExperiment
from MSChromatogram cimport MSChromatogram
from MzTabMetaData cimport MzTabMetaData
from String cimport String
from Types cimport UInt


cdef extern from "<OpenMS/QC/TIC.h>" namespace "OpenMS":

    cdef cppclass TIC(QCBase):

        # nested struct
        cdef cppclass Result:
            vector[UInt] intensities
            vector[float] relative_intensities
            vector[float] retention_times
            UInt area
            UInt fall
            UInt jump

        TIC() except +

        Result compute(const MSExperiment& exp,
                       float bin_size = 0,
                       UInt ms_level = 1) except +

        const String& getName() const except +

        const vector[MSChromatogram]& getResults() const except +

        QCBase.Status requirements() const except +

        void addMetaDataMetricsToMzTab(MzTabMetaData& meta,
                                       vector[Result]& tics) except +