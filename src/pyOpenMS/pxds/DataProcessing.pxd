from Types cimport *
from SpectrumSettings cimport *
from Peak1D cimport *
from String cimport *
from Software cimport *
from DateTime cimport *
from MetaInfoInterface cimport *

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

    ctypedef shared_ptr[DataProcessing] DataProcessingPtr


cdef extern from "<OpenMS/METADATA/DataProcessing.h>" namespace "OpenMS::DataProcessing":

    cdef enum ProcessingAction:

        DATA_PROCESSING, CHARGE_DECONVOLUTION, DEISOTOPING, SMOOTHING,
        CHARGE_CALCULATION, PRECURSOR_RECALCULATION, BASELINE_REDUCTION,
        PEAK_PICKING, ALIGNMENT, CALIBRATION, NORMALIZATION, FILTERING,
        QUANTITATION, FEATURE_GROUPING, IDENTIFICATION_MAPPING,
        FORMAT_CONVERSION, CONVERSION_MZDATA, CONVERSION_MZML,
        CONVERSION_MZXML, CONVERSION_DTA, SIZE_OF_PROCESSINGACTION

