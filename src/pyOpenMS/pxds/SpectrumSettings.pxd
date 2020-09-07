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

        SpectrumSettings()    nogil except +
        void unify(SpectrumSettings)    nogil except +
        int  getType()    nogil except +
        void setType(SpectrumType)    nogil except +
        String getNativeID()    nogil except +
        void setNativeID(String)    nogil except +
        String getComment()    nogil except +
        void setComment(String)    nogil except +

        InstrumentSettings getInstrumentSettings()    nogil except +
        void setInstrumentSettings(InstrumentSettings)    nogil except +

        AcquisitionInfo getAcquisitionInfo()    nogil except +
        void setAcquisitionInfo(AcquisitionInfo)    nogil except +

        SourceFile getSourceFile()    nogil except +
        void setSourceFile(SourceFile)    nogil except +

        libcpp_vector[Precursor] getPrecursors() nogil except +
        void setPrecursors(libcpp_vector[Precursor])   nogil except +

        libcpp_vector[Product] getProducts() nogil except +
        void setProducts(libcpp_vector[Product])   nogil except +

        libcpp_vector[PeptideIdentification] getPeptideIdentifications() nogil except +
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification])   nogil except +

        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() nogil except +
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ])   nogil except +

cdef extern from "<OpenMS/METADATA/SpectrumSettings.h>" namespace "OpenMS::SpectrumSettings":

    cdef enum SpectrumType:
        # wrap-attach:
        #     SpectrumSettings
        UNKNOWN,  PEAKS, RAWDATA, SIZE_OF_SPECTRUMTYPE

