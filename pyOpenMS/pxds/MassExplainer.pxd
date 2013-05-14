from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from Adduct cimport *
from Compomer cimport *

ctypedef libcpp_vector[Adduct] AdductsType

cdef extern from "<OpenMS/DATASTRUCTURES/MassExplainer.h>" namespace "OpenMS":
    
    cdef cppclass MassExplainer "OpenMS::MassExplainer":
        MassExplainer() nogil except +
        MassExplainer(MassExplainer) nogil except + #wrap-ignore
        ## MassExplainer(libcpp_vector[Adduct] adduct_base) nogil except +
        ## MassExplainer(Int q_min, Int q_max, Int max_span, DoubleReal thresh_logp) nogil except +
        ## MassExplainer(libcpp_vector[Adduct] adduct_base, Int q_min, Int q_max, Int max_span, DoubleReal thresh_logp, Size max_neutrals) nogil except +
        ## ## MassExplainer operator=(MassExplainer &rhs) nogil except +
        void setAdductBase(libcpp_vector[Adduct] adduct_base) nogil except +
        libcpp_vector[Adduct] getAdductBase() nogil except +
        Compomer getCompomerById(Size id) nogil except +
        void compute() nogil except +

