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
        #   MetaInfoInterface
        #
        # wrap-doc:
        #  Description of the chromatogram settings, provides meta-information
        #  about a single chromatogram.

        ChromatogramSettings() except + nogil 
        ChromatogramSettings(ChromatogramSettings &) except + nogil 

        Product getProduct() except + nogil  # wrap-doc:Returns the product ion
        void setProduct(Product p) except + nogil  # wrap-doc:Sets the product ion

        String getNativeID() except + nogil  # wrap-doc:Returns the native identifier for the spectrum, used by the acquisition software.
        void setNativeID(String native_id) except + nogil  # wrap-doc:Sets the native identifier for the spectrum, used by the acquisition software.

        # returns the free-text comment
        String getComment() except + nogil  # wrap-doc:Returns the free-text comment
        # sets the free-text comment
        void setComment(String comment) except + nogil  # wrap-doc:Sets the free-text comment

        # returns a mutable reference to the instrument settings of the current spectrum
        InstrumentSettings getInstrumentSettings() except + nogil  # wrap-doc:Returns the instrument settings of the current spectrum
        # sets the instrument settings of the current spectrum
        void setInstrumentSettings(InstrumentSettings instrument_settings) except + nogil  # wrap-doc:Sets the instrument settings of the current spectrum

        # returns a mutable reference to the acquisition info
        AcquisitionInfo getAcquisitionInfo() except + nogil  # wrap-doc:Returns the acquisition info
        # sets the acquisition info
        void setAcquisitionInfo(AcquisitionInfo acquisition_info) except + nogil  # wrap-doc:Sets the acquisition info

        # returns a mutable reference to the source file
        SourceFile getSourceFile() except + nogil  # wrap-doc:Returns the source file
        # sets the source file
        void setSourceFile(SourceFile source_file) except + nogil  # wrap-doc:Sets the source file

        # returns a mutable reference to the precursors
        Precursor getPrecursor() except + nogil  # wrap-doc:Returns the precursors
        # sets the precursors
        void setPrecursor(Precursor precursor) except + nogil  # wrap-doc:Sets the precursors

        # returns a mutable reference to the description of the applied processing
        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() except + nogil  # wrap-doc:Returns the description of the applied processing
        # sets the description of the applied processing
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ])   except + nogil  # wrap-doc:Sets the description of the applied processing

        # sets the chromatogram type
        void setChromatogramType(ChromatogramType type) except + nogil  # wrap-doc:Sets the chromatogram type
        ChromatogramType getChromatogramType() except + nogil  # wrap-doc:Get the chromatogram type


cdef extern from "<OpenMS/METADATA/ChromatogramSettings.h>" namespace "OpenMS::ChromatogramSettings":

    cdef enum ChromatogramType:
        # wrap-attach:
        #    ChromatogramSettings

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
