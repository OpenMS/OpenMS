from Types cimport *
from Peak1D cimport *
from String cimport *
from Software cimport *
from DateTime cimport *
from MetaInfoInterface cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/METADATA/DataProcessing.h>" namespace "OpenMS":

    cdef cppclass DataProcessing(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        DataProcessing()  except + nogil 
        DataProcessing(DataProcessing &) except + nogil 

        void setProcessingActions(libcpp_set[ProcessingAction]) except + nogil 
        libcpp_set[ProcessingAction] getProcessingActions() except + nogil 

        Software getSoftware() except + nogil 
        void setSoftware(Software s) except + nogil 

        DateTime getCompletionTime()  except + nogil 
        void setCompletionTime(DateTime t) except + nogil 

        @staticmethod
        libcpp_vector[String] getAllNamesOfProcessingAction() except + nogil  # wrap-doc:Returns all processing action names known to OpenMS

        @staticmethod
        String processingActionToString(ProcessingAction action) except + nogil  # wrap-doc:Convert a ProcessingAction enum to String. Throws Exception::InvalidValue if action is SIZE_OF_PROCESSINGACTION

        @staticmethod
        ProcessingAction toProcessingAction(const String& name) except + nogil  # wrap-doc:Convert a string to ProcessingAction enum. Throws Exception::InvalidValue if name is not found

    ctypedef shared_ptr[DataProcessing] DataProcessingPtr


cdef extern from "<OpenMS/METADATA/DataProcessing.h>" namespace "OpenMS::DataProcessing":

    cdef enum ProcessingAction:
        # wrap-attach:
        #   DataProcessing
        DATA_PROCESSING,                #< General data processing (if no other term applies)
        CHARGE_DECONVOLUTION,           #< Charge deconvolution
        DEISOTOPING,                    #< Deisotoping
        SMOOTHING,                      #< Smoothing of the signal to reduce noise
        CHARGE_CALCULATION,             #< Determination of the peak charge
        PRECURSOR_RECALCULATION,        #< Recalculation of precursor m/z
        BASELINE_REDUCTION,             #< Baseline reduction
        PEAK_PICKING,                   #< Peak picking (conversion from raw to peak data)
        ALIGNMENT,                      #< Retention time alignment of different maps
        CALIBRATION,                    #< Calibration of m/z positions
        NORMALIZATION,                  #< Normalization of intensity values
        FILTERING,                      #< Data filtering or extraction
        QUANTITATION,                   #< Quantitation
        FEATURE_GROUPING,               #< Feature grouping
        IDENTIFICATION_MAPPING,         #< Identification mapping
        FORMAT_CONVERSION,              #< General file format conversion (if no other term applies)
        CONVERSION_MZDATA,              #< Conversion to mzData format
        CONVERSION_MZML,                #< Conversion to mzML format
        CONVERSION_MZXML,               #< Conversion to mzXML format
        CONVERSION_DTA,                 #< Conversion to DTA format
        IDENTIFICATION,                 #< Identification
        ION_MOBILITY_BINNING,           #< Ion mobility binning (merging of spectra with similar IM values)
        SIZE_OF_PROCESSINGACTION
