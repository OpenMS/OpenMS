from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from String cimport *
from QCBase cimport *
from MSExperiment cimport *
from MSChromatogram cimport *

cdef extern from "<OpenMS/QC/TIC.h>" namespace "OpenMS":
    
    cdef cppclass TIC(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  Total Ion Count (TIC) as a QC metric
        #  
        #  Simple class to calculate the TIC of an MSExperiment.
        #  Allows for multiple usage, because each calculated TIC is
        #  stored internally. Those results can then be returned using
        #  getResults().

        TIC() except + nogil 
        TIC(TIC &) except + nogil 

        Result compute(MSExperiment & exp, float bin_size, UInt ms_level) except + nogil 
            # wrap-doc:
            #  Compute Total Ion Count and applies the resampling algorithm, 
            #  if a bin size in RT seconds greater than 0 is given.
            #  All MS1 TICs within a bin are summed up.
            #  
            #  :param exp: Peak map to compute the MS1 tick from
            #  :param bin_size: RT bin size in seconds (default 0)
            #  :param ms_level: MS level of spectra for calculation (default 1)
            #  :return: Result struct with computed QC metrics

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric

        libcpp_vector[MSChromatogram] getResults() except + nogil  # wrap-doc:Returns results

cdef extern from "<OpenMS/QC/TIC.h>" namespace "OpenMS::TIC":
    
    cdef cppclass Result "OpenMS::TIC::Result":
        # wrap-doc:
        #  Structure for storing TIC computation results
        
        Result() except + nogil 
        Result(Result &) except + nogil 
        
        libcpp_vector[UInt] intensities  # wrap-doc:TIC intensities
        libcpp_vector[float] relative_intensities  # wrap-doc:Relative TIC intensities
        libcpp_vector[float] retention_times  # wrap-doc:TIC retention times in seconds
        UInt area  # wrap-doc:Area under TIC
        UInt fall  # wrap-doc:MS1 signal fall (10x) count
        UInt jump  # wrap-doc:MS1 signal jump (10x) count
        
        bool operator==(Result & rhs) except + nogil 
