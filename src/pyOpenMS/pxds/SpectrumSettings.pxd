from Types cimport *
from String cimport *
from Peak1D cimport *
from InstrumentSettings cimport *
from SourceFile cimport *
from PeptideIdentification cimport *
from Precursor cimport *
from DataProcessing cimport *
from Product cimport *
from AcquisitionInfo cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/SpectrumSettings.h>" namespace "OpenMS":

    cdef cppclass SpectrumSettings(MetaInfoInterface):
        # wrap-inherits:
        #  MetaInfoInterface

        SpectrumSettings() except + nogil 
        SpectrumSettings(SpectrumSettings &) except + nogil 

        void unify(SpectrumSettings) except + nogil 
        int  getType() except + nogil  # wrap-doc:Returns the spectrum type (centroided (PEAKS) or profile data (RAW))
        void setType(SpectrumType) except + nogil  # wrap-doc:Sets the spectrum type
        String getNativeID() except + nogil  # wrap-doc:Returns the native identifier for the spectrum, used by the acquisition software
        void setNativeID(String) except + nogil  # wrap-doc:Sets the native identifier for the spectrum, used by the acquisition software
        String getComment() except + nogil  # wrap-doc:Returns the free-text comment
        void setComment(String) except + nogil  # wrap-doc:Sets the free-text comment

        InstrumentSettings getInstrumentSettings() except + nogil  # wrap-doc:Returns a const reference to the instrument settings of the current spectrum
        void setInstrumentSettings(InstrumentSettings) except + nogil  # wrap-doc:Sets the instrument settings of the current spectrum

        AcquisitionInfo getAcquisitionInfo() except + nogil  # wrap-doc:Returns a const reference to the acquisition info
        void setAcquisitionInfo(AcquisitionInfo) except + nogil  # wrap-doc:Sets the acquisition info

        SourceFile getSourceFile() except + nogil  # wrap-doc:Returns a const reference to the source file
        void setSourceFile(SourceFile) except + nogil  # wrap-doc:Sets the source file

        libcpp_vector[Precursor] getPrecursors() except + nogil  # wrap-doc:Returns a const reference to the precursors
        void setPrecursors(libcpp_vector[Precursor]) except + nogil  # wrap-doc:Sets the precursors

        libcpp_vector[Product] getProducts() except + nogil  # wrap-doc:Returns a const reference to the products
        void setProducts(libcpp_vector[Product]) except + nogil  # wrap-doc:Sets the products

        libcpp_vector[PeptideIdentification] getPeptideIdentifications() except + nogil  # wrap-doc:Returns a const reference to the PeptideIdentification vector
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification]) except + nogil  # wrap-doc:Sets the PeptideIdentification vector

        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() except + nogil 
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ]) except + nogil 

cdef extern from "<OpenMS/METADATA/SpectrumSettings.h>" namespace "OpenMS::SpectrumSettings":

    cdef enum SpectrumType:
        # wrap-attach:
        #    SpectrumSettings
        UNKNOWN, CENTROID, PROFILE, SIZE_OF_SPECTRUMTYPE
