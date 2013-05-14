from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from Adduct cimport *
from StringList cimport *
from Map cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Compomer.h>" namespace "OpenMS":

    cdef cppclass Compomer:
  
        Compomer() nogil except +
        Compomer(Compomer) nogil except + 
  
        void add(Adduct & a, UInt side) nogil except +

        bool isConflicting(Compomer & cmp, UInt side_this, UInt side_other) nogil except +

        # /// set an Id which allows unique identification of a compomer
        void setID(Size id) nogil except +

        # /// return Id which allows unique identification of this compomer
        Size getID() nogil except +

        # /// left and right adducts of this compomer
        # TODO OpenMS Map type
        libcpp_vector[Map[String, Adduct] ] getComponent() nogil except + # wrap-ignore
    
        # /// net charge of compomer (i.e. difference between left and right side of compomer)
        Int getNetCharge() nogil except +

        # /// mass of all contained adducts
        DoubleReal getMass() nogil except +

        # /// summed positive charges of contained adducts
        Int getPositiveCharges() nogil except +

        # /// summed negative charges of contained adducts
        Int getNegativeCharges() nogil except +

        # /// return log probability
        DoubleReal getLogP() nogil except +

        # /// return log probability
        DoubleReal getRTShift() nogil except +

        # /// get adducts with their abundance as compact string for both sides
        String getAdductsAsString() nogil except +

        # /// get adducts with their abundance as compact string (amounts are absolute unless side=BOTH)
        # /// @param side Use LEFT for left, RIGHT for right
        String getAdductsAsString(UInt side) nogil except +

        # /// check if Compomer only contains a single adduct on side @p side
        bool isSingleAdduct(Adduct & a, UInt side) nogil except +

        Compomer removeAdduct(Adduct & a) nogil except +

        Compomer removeAdduct(Adduct & a, UInt side) nogil except +

        StringList getLabels(UInt side) nogil except +

        # void add(CompomerSide & add_side, UInt side) nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/Compomer.h>" namespace "OpenMS::Compomer":
    # side of compomer (LEFT ^ substract; RIGHT ^ add)
    cdef enum SIDE: LEFT, RIGHT, BOTH
