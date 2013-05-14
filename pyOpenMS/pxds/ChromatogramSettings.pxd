from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from InstrumentSettings cimport *
from Precursor cimport *
from Peak1D cimport *
from SourceFile cimport *
from PeptideIdentification cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from DataProcessing cimport *

from Product cimport *
from AcquisitionInfo cimport *
from InstrumentSettings cimport *
from SourceFile cimport *

cdef extern from "<OpenMS/METADATA/ChromatogramSettings.h>" namespace "OpenMS":

    cdef cppclass ChromatogramSettings:

        ChromatogramSettings() nogil except +
        ChromatogramSettings(ChromatogramSettings) nogil except +

        Product getProduct() nogil except +
        void setProduct(Product p) nogil except +

        # returns the native identifier for the spectrum, used by the acquisition software.
        String getNativeID() nogil except +
        # sets the native identifier for the spectrum, used by the acquisition software.
        void setNativeID(String native_id) nogil except +

        # returns the free-text comment
        String getComment() nogil except +
        # sets the free-text comment
        void setComment(String comment) nogil except +

        # returns a mutable reference to the instrument settings of the current spectrum
        InstrumentSettings getInstrumentSettings() nogil except +
        # sets the instrument settings of the current spectrum
        void setInstrumentSettings(InstrumentSettings instrument_settings) nogil except +

        # returns a mutable reference to the acquisition info
        AcquisitionInfo getAcquisitionInfo() nogil except +
        # sets the acquisition info
        void setAcquisitionInfo(AcquisitionInfo acquisition_info) nogil except +

        # returns a mutable reference to the source file
        SourceFile getSourceFile() nogil except +
        # sets the source file
        void setSourceFile(SourceFile source_file) nogil except +

        # returns a mutable reference to the precursors
        Precursor getPrecursor() nogil except +
        # sets the precursors
        void setPrecursor(Precursor precursor) nogil except +

        # returns a reference to the description of the applied processing
        libcpp_vector[DataProcessing] getDataProcessing() nogil except +
        # sets the description of the applied processing
        void setDataProcessing(libcpp_vector[DataProcessing]) nogil except +

        # # returns the chromatogram type, e.g. a SRM chromatogram
        # ChromatogramType getChromatogramType() nogil except +
        # # sets the chromatogram type
        # void setChromatogramType(ChromatogramType type) nogil except +

