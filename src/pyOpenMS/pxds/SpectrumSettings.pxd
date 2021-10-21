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
        #   MetaInfoInterface

        SpectrumSettings() nogil except +
        SpectrumSettings(SpectrumSettings &) nogil except +

        void unify(SpectrumSettings) nogil except +
        int  getType() nogil except + # wrap-doc:Returns the spectrum type (centroided (PEAKS) or profile data (RAW))
        void setType(SpectrumType) nogil except + # wrap-doc:Sets the spectrum type
        String getNativeID() nogil except + # wrap-doc:Returns the native identifier for the spectrum, used by the acquisition software
        void setNativeID(String) nogil except + # wrap-doc:Sets the native identifier for the spectrum, used by the acquisition software
        String getComment() nogil except + # wrap-doc:Returns the free-text comment
        void setComment(String) nogil except + # wrap-doc:Sets the free-text comment

        InstrumentSettings getInstrumentSettings() nogil except + # wrap-doc:Returns a const reference to the instrument settings of the current spectrum
        void setInstrumentSettings(InstrumentSettings) nogil except + # wrap-doc:Sets the instrument settings of the current spectrum

        AcquisitionInfo getAcquisitionInfo() nogil except + # wrap-doc:Returns a const reference to the acquisition info
        void setAcquisitionInfo(AcquisitionInfo) nogil except + # wrap-doc:Sets the acquisition info

        SourceFile getSourceFile() nogil except + # wrap-doc:Returns a const reference to the source file
        void setSourceFile(SourceFile) nogil except + # wrap-doc:Sets the source file

        libcpp_vector[Precursor] getPrecursors() nogil except + # wrap-doc:Returns a const reference to the precursors
        void setPrecursors(libcpp_vector[Precursor]) nogil except + # wrap-doc:Sets the precursors

        libcpp_vector[Product] getProducts() nogil except + # wrap-doc:Returns a const reference to the products
        void setProducts(libcpp_vector[Product]) nogil except + # wrap-doc:Sets the products

        libcpp_vector[PeptideIdentification] getPeptideIdentifications() nogil except + # wrap-doc:Returns a const reference to the PeptideIdentification vector
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification]) nogil except + # wrap-doc:Sets the PeptideIdentification vector

        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() nogil except +
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ]) nogil except +

cdef extern from "<OpenMS/METADATA/SpectrumSettings.h>" namespace "OpenMS::SpectrumSettings":

    cdef enum SpectrumType:
        # wrap-attach:
        #     SpectrumSettings
        UNKNOWN,  PEAKS, RAWDATA, SIZE_OF_SPECTRUMTYPE
