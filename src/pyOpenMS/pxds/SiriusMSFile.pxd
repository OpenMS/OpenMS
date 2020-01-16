from Types cimport *
from DefaultParamHandler cimport *
from BaseFeature cimport *
from ProgressLogger cimport *
from String cimport *
from StringList cimport *
from MSExperiment cimport *
from Map cimport *
from FeatureMapping cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusMSConverter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMSFile "OpenMS::SiriusMSFile":

        SiriusMSFile() nogil except +
        SiriusMSFile(SiriusMSFile) nogil except + #wrap-ignore

        void store(MSExperiment& spectra,
                   String& msfile,
                   FeatureMapping_FeatureToMs2Indices& feature_ms2_spectra_map,
                   bool& feature_only,
                   int& isotope_pattern_iterations,
                   bool no_mt_info,
                   libcpp_vector[ SiriusMSFile_CompoundInfo ] v_cmpinfo) nogil except +

    cdef cppclass SiriusMSFile_CompoundInfo "OpenMS::SiriusMSFile::CompoundInfo":

        SiriusMSFile_CompoundInfo() nogil except +
        SiriusMSFile_CompoundInfo(SiriusMSFile_CompoundInfo) nogil except + #wrap-ignore

    cdef cppclass SiriusMSFile_AccessionInfo "OpenMS::SiriusMSFile::AccessionInfo":

        SiriusMSFile_AccessionInfo() nogil except +
        SiriusMSFile_AccessionInfo(SiriusMSFile_AccessionInfo) nogil except + #wrap-ignore