from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Adduct.h>" namespace "OpenMS":

    cdef cppclass Adduct:
  
        Adduct() nogil except +
        Adduct(Adduct) nogil except + 
  
        Adduct(Int charge) nogil except +
        Adduct(Int charge, Int amount, double singleMass, String formula, double log_prob, double rt_shift, String label) nogil except +
  
        Int getCharge() nogil except +
  
        void setCharge(Int charge) nogil except +
  
        Int getAmount() nogil except +
        void setAmount(Int amount) nogil except +
  
        double getSingleMass() nogil except +
        void setSingleMass(double singleMass) nogil except +
  
        double getLogProb() nogil except +
        void setLogProb(double log_prob) nogil except +
  
        String getFormula() nogil except +
        void setFormula(String formula) nogil except +
  
        double getRTShift() nogil except +
        String getLabel() nogil except +
