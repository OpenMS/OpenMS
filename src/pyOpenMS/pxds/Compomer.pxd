from libcpp cimport bool
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from Adduct cimport *
from StringList cimport *
from Map cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Compomer.h>" namespace "OpenMS":

    cdef cppclass Compomer:
  
        Compomer() nogil except +
        Compomer(Compomer &) nogil except +
  
        void add(Adduct & a, UInt side) nogil except +

        bool isConflicting(Compomer & cmp, UInt side_this, UInt side_other) nogil except +

        # /// set an Id which allows unique identification of a compomer
        void setID(Size id) nogil except + # wrap-doc:Sets an Id which allows unique identification of a compomer

        # /// return Id which allows unique identification of this compomer
        Size getID() nogil except + # wrap-doc:Returns Id which allows unique identification of this compomer

        # /// left and right adducts of this compomer
        # TODO OpenMS Map type
        libcpp_vector[Map[String, Adduct] ] getComponent() nogil except + # wrap-ignore
    
        # /// net charge of compomer (i.e. difference between left and right side of compomer)
        Int getNetCharge() nogil except + # wrap-doc:Net charge of compomer (i.e. difference between left and right side of compomer)

        # /// mass of all contained adducts
        double getMass() nogil except + # wrap-doc:Mass of all contained adducts

        # /// summed positive charges of contained adducts
        Int getPositiveCharges() nogil except + # wrap-doc:Summed positive charges of contained adducts

        # /// summed negative charges of contained adducts
        Int getNegativeCharges() nogil except + # wrap-doc:Summed negative charges of contained adducts

        # /// return log probability
        double getLogP() nogil except + # wrap-doc:Returns the log probability

        # /// return log probability
        double getRTShift() nogil except + # wrap-doc:Returns the log probability

        # /// get adducts with their abundance as compact string for both sides
        String getAdductsAsString() nogil except + # wrap-doc:Get adducts with their abundance as compact string for both sides

        # /// get adducts with their abundance as compact string (amounts are absolute unless side=BOTH)
        # /// @param side Use LEFT for left, RIGHT for right
        String getAdductsAsString(UInt side) nogil except + # wrap-doc:Get adducts with their abundance as compact string (amounts are absolute unless side=BOTH)

        # /// check if Compomer only contains a single adduct on side @p side
        bool isSingleAdduct(Adduct & a, UInt side) nogil except + # wrap-doc:Check if Compomer only contains a single adduct on side @p side

        Compomer removeAdduct(Adduct & a) nogil except + # wrap-doc:Remove ALL instances of the given adduct

        Compomer removeAdduct(Adduct & a, UInt side) nogil except +

        StringList getLabels(UInt side) nogil except + # wrap-doc:Returns the adduct labels from parameter(side) given. (LEFT or RIGHT)

        # void add(CompomerSide & add_side, UInt side) nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/Compomer.h>" namespace "OpenMS::Compomer":
    # side of compomer (LEFT ^ substract; RIGHT ^ add)
    cdef enum SIDE: LEFT, RIGHT, BOTH
