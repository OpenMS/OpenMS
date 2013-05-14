from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Adduct.h>" namespace "OpenMS":

    cdef cppclass Adduct:
  
        Adduct() nogil except +
        Adduct(Adduct) nogil except + 
  
        Adduct(Int charge) nogil except +
        Adduct(Int charge, Int amount, DoubleReal singleMass, String formula, DoubleReal log_prob, DoubleReal rt_shift, String label) nogil except +
  
        Int getCharge() nogil except +
  
        void setCharge(Int charge) nogil except +
  
        Int getAmount() nogil except +
        void setAmount(Int amount) nogil except +
  
        DoubleReal getSingleMass() nogil except +
        void setSingleMass(DoubleReal singleMass) nogil except +
  
        DoubleReal getLogProb() nogil except +
        void setLogProb(DoubleReal log_prob) nogil except +
  
        String getFormula() nogil except +
        void setFormula(String formula) nogil except +
  
        DoubleReal getRTShift() nogil except +
        String getLabel() nogil except +
