from Types cimport *
from String cimport *
from MSExperiment cimport *
from FeatureMapping cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusMSConverter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMSFile:
        SiriusMSFile() nogil except +
        SiriusMSFile(SiriusMSFile) nogil except + #wrap-ignore

    cdef cppclass SiriusMSFile_CompoundInfo "OpenMS::SiriusMSFile::CompoundInfo":
        SiriusMSFile_CompoundInfo() nogil except +
        SiriusMSFile_CompoundInfo(SiriusMSFile_CompoundInfo) nogil except + #wrap-ignore

    cdef cppclass SiriusMSFile_AccessionInfo "OpenMS::SiriusMSFile::AccessionInfo":
        SiriusMSFile_AccessionInfo() nogil except +
        SiriusMSFile_AccessionInfo(SiriusMSFile_AccessionInfo) nogil except + #wrap-ignore

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusMSConverter.h>" namespace "OpenMS::SiriusMSFile":

        # wrap static method:
        void store(MSExperiment spectra,
                   String msfile,
                   FeatureMapping_FeatureToMs2Indices feature_mapping,
                   bool feature_only,
                   int isotope_pattern_iterations,
                   bool no_mt_info,
                   libcpp_vector[ SiriusMSFile_CompoundInfo ] v_cmpinfo) nogil except + # wrap-attach:SiriusMSFile
