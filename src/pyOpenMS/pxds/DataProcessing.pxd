from Types cimport *
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

        DATA_PROCESSING,                # General data processing (if no other term applies)
        CHARGE_DECONVOLUTION,           # Charge deconvolution
        DEISOTOPING,                    # Deisotoping
        SMOOTHING,                      # Smoothing of the signal to reduce noise
        CHARGE_CALCULATION,             # Determination of the peak charge
        PRECURSOR_RECALCULATION,        # Recalculation of precursor m/z
        BASELINE_REDUCTION,             # Baseline reduction
        PEAK_PICKING,                   # Peak picking (conversion from raw to peak data)
        ALIGNMENT,                      # Retention time alignment of different maps
        CALIBRATION,                    # Calibration of m/z positions
        NORMALIZATION,                  # Normalization of intensity values
        FILTERING,                      # Data filtering or extraction
        QUANTITATION,                   # Quantitation
        FEATURE_GROUPING,               # %Feature grouping
        IDENTIFICATION_MAPPING,         # %Identification mapping
        FORMAT_CONVERSION,              # General file format conversion (if no other term applies)
        CONVERSION_MZDATA,              # Conversion to mzData format
        CONVERSION_MZML,                # Conversion to mzML format
        CONVERSION_MZXML,               # Conversion to mzXML format
        CONVERSION_DTA,                 # Conversion to DTA format
        SIZE_OF_PROCESSINGACTION

