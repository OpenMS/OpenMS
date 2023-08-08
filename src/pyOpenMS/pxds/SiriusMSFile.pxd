from Types cimport *
from String cimport *
from FeatureMapping cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusMSConverter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMSFile:
        SiriusMSFile() except + nogil 
        SiriusMSFile(SiriusMSFile &) except + nogil  # compiler

    cdef cppclass SiriusMSFile_CompoundInfo "OpenMS::SiriusMSFile::CompoundInfo":
        SiriusMSFile_CompoundInfo() except + nogil 
        SiriusMSFile_CompoundInfo(SiriusMSFile_CompoundInfo &) except + nogil  # compiler

    cdef cppclass SiriusMSFile_AccessionInfo "OpenMS::SiriusMSFile::AccessionInfo":
        SiriusMSFile_AccessionInfo() except + nogil 
        SiriusMSFile_AccessionInfo(SiriusMSFile_AccessionInfo &) except + nogil  # compiler

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusMSConverter.h>" namespace "OpenMS::SiriusMSFile":

        # wrap static method:
        void store(MSExperiment& spectra,
                   const String& msfile,
                   FeatureMapping_FeatureToMs2Indices& feature_mapping,
                   bool feature_only,
                   int isotope_pattern_iterations,
                   bool no_mt_info,
                   libcpp_vector[ SiriusMSFile_CompoundInfo ]& v_cmpinfo) except + nogil  # wrap-attach:SiriusMSFile
