from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from Adduct cimport *
from Compomer cimport *

ctypedef libcpp_vector[Adduct] AdductsType

cdef extern from "<OpenMS/DATASTRUCTURES/MassExplainer.h>" namespace "OpenMS":
    
    cdef cppclass MassExplainer "OpenMS::MassExplainer":
        MassExplainer() nogil except + # wrap-doc:Computes empirical formulas for given mass differences using a set of allowed elements
        MassExplainer(MassExplainer &) nogil except + # compiler

        MassExplainer(libcpp_vector[Adduct] adduct_base) nogil except +
        MassExplainer(Int q_min, Int q_max, Int max_span, double thresh_logp) nogil except +
        ## MassExplainer(libcpp_vector[Adduct] adduct_base, Int q_min, Int q_max, Int max_span, double thresh_logp, Size max_neutrals) nogil except +
        ## MassExplainer operator=(MassExplainer &rhs) nogil except +
        void setAdductBase(libcpp_vector[Adduct] adduct_base) nogil except + # wrap-doc:Sets the set of possible adducts
        libcpp_vector[Adduct] getAdductBase() nogil except + # wrap-doc:Returns the set of adducts
        Compomer getCompomerById(Size id) nogil except + # wrap-doc:Returns a compomer by its Id (useful after a query() )
        void compute() nogil except + # wrap-doc:Fill map with possible mass-differences along with their explanation
