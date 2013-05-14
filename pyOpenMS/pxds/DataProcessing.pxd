from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp.string cimport string as libcpp_string
from SpectrumSettings cimport *
from Peak1D cimport *
from String cimport *
from Software cimport *
from DateTime cimport *

cdef extern from "<OpenMS/METADATA/DataProcessing.h>" namespace "OpenMS":

    cdef cppclass DataProcessing(MetaInfoInterface):
        # wrap-inherits:
        #     MetaInfoInterface

        DataProcessing()  nogil except +
        DataProcessing(DataProcessing) nogil except + # wrap-ignore

        void setProcessingActions(libcpp_set[ProcessingAction]) nogil except +
        libcpp_set[ProcessingAction] getProcessingActions() nogil except +

        Software getSoftware() nogil except +
        void setSoftware(Software s) nogil except +

        DateTime getCompletionTime()  nogil except +
        void setCompletionTime(DateTime t) nogil except +

        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:

        void getKeys(libcpp_vector[String] & keys)
        void getKeys(libcpp_vector[unsigned int] & keys)
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +


cdef extern from "<OpenMS/METADATA/DataProcessing.h>" namespace "OpenMS::DataProcessing":

    cdef enum ProcessingAction:

        DATA_PROCESSING, CHARGE_DECONVOLUTION, DEISOTOPING, SMOOTHING,
        CHARGE_CALCULATION, PRECURSOR_RECALCULATION, BASELINE_REDUCTION,
        PEAK_PICKING, ALIGNMENT, CALIBRATION, NORMALIZATION, FILTERING,
        QUANTITATION, FEATURE_GROUPING, IDENTIFICATION_MAPPING,
        FORMAT_CONVERSION, CONVERSION_MZDATA, CONVERSION_MZML,
        CONVERSION_MZXML, CONVERSION_DTA, SIZE_OF_PROCESSINGACTION
