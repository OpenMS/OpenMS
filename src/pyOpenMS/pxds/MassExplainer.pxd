from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from Adduct cimport *
from Compomer cimport *

ctypedef libcpp_vector[Adduct] AdductsType

cdef extern from "<OpenMS/DATASTRUCTURES/MassExplainer.h>" namespace "OpenMS":
    
    cdef cppclass MassExplainer "OpenMS::MassExplainer":
        MassExplainer() except + nogil  # wrap-doc:Computes empirical formulas for given mass differences using a set of allowed elements
        MassExplainer(MassExplainer &) except + nogil  # compiler

        MassExplainer(libcpp_vector[Adduct] adduct_base) except + nogil 
        MassExplainer(Int q_min, Int q_max, Int max_span, double thresh_logp) except + nogil 
        ## MassExplainer(libcpp_vector[Adduct] adduct_base, Int q_min, Int q_max, Int max_span, double thresh_logp, Size max_neutrals) except + nogil 
        ## MassExplainer operator=(MassExplainer &rhs) except + nogil 
        void setAdductBase(libcpp_vector[Adduct] adduct_base) except + nogil  # wrap-doc:Sets the set of possible adducts
        libcpp_vector[Adduct] getAdductBase() except + nogil  # wrap-doc:Returns the set of adducts
        Compomer getCompomerById(Size id) except + nogil  # wrap-doc:Returns a compomer by its Id (useful after a query() )
        void compute() except + nogil  # wrap-doc:Fill map with possible mass-differences along with their explanation
