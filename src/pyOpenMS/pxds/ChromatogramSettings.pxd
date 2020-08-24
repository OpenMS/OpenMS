from Types cimport *
from String cimport *
from InstrumentSettings cimport *
from SourceFile cimport *
from PeptideIdentification cimport *
from DataProcessing cimport *
from Product cimport *
from AcquisitionInfo cimport *
from Precursor cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/ChromatogramSettings.h>" namespace "OpenMS":

    cdef cppclass ChromatogramSettings(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface
        #
        # wrap-doc:
        #   Description of the chromatogram settings, provides meta-information
        #   about a single chromatogram.

        ChromatogramSettings() nogil except +
        ChromatogramSettings(ChromatogramSettings) nogil except +

        Product getProduct() nogil except + # wrap-doc:Access to the product ion
        void setProduct(Product p) nogil except + # wrap-doc:Set the product ion

        String getNativeID() nogil except + # wrap-doc:returns the native identifier for the spectrum, used by the acquisition software.
        void setNativeID(String native_id) nogil except + # wrap-doc:sets the native identifier for the spectrum, used by the acquisition software.

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

        # returns a mutable reference to the description of the applied processing
        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() nogil except +
        # sets the description of the applied processing
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ])   nogil except +

        # sets the chromatogram type
        void setChromatogramType(ChromatogramType type) nogil except +
        ChromatogramType getChromatogramType() nogil except +


cdef extern from "<OpenMS/METADATA/ChromatogramSettings.h>" namespace "OpenMS::ChromatogramSettings":

    cdef enum ChromatogramType:
        # wrap-attach:
        #     ChromatogramSettings

        MASS_CHROMATOGRAM,
        TOTAL_ION_CURRENT_CHROMATOGRAM,
        SELECTED_ION_CURRENT_CHROMATOGRAM,
        BASEPEAK_CHROMATOGRAM,
        SELECTED_ION_MONITORING_CHROMATOGRAM,
        SELECTED_REACTION_MONITORING_CHROMATOGRAM,
        ELECTROMAGNETIC_RADIATION_CHROMATOGRAM,
        ABSORPTION_CHROMATOGRAM,
        EMISSION_CHROMATOGRAM,
        SIZE_OF_CHROMATOGRAM_TYPE

