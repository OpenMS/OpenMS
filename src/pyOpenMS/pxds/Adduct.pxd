from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Adduct.h>" namespace "OpenMS":

    cdef cppclass Adduct:
  
        Adduct() except + nogil 
        Adduct(Adduct &) except + nogil  # compiler
  
        Adduct(Int charge) except + nogil 
        Adduct(Int charge, Int amount, double singleMass, String formula, double log_prob, double rt_shift, String label) except + nogil 
  
        Int getCharge() except + nogil 
  
        void setCharge(Int charge) except + nogil 
  
        Int getAmount() except + nogil 
        void setAmount(Int amount) except + nogil 
  
        double getSingleMass() except + nogil 
        void setSingleMass(double singleMass) except + nogil 
  
        double getLogProb() except + nogil 
        void setLogProb(double log_prob) except + nogil 
  
        String getFormula() except + nogil 
        void setFormula(String formula) except + nogil 
  
        double getRTShift() except + nogil 
        String getLabel() except + nogil 
