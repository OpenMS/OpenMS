from Types cimport *
from Adduct cimport *
from StringList cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Compomer.h>" namespace "OpenMS":

    cdef cppclass Compomer:
  
        Compomer() except + nogil 
        Compomer(Compomer &) except + nogil 
  
        void add(Adduct & a, UInt side) except + nogil 

        bool isConflicting(Compomer & cmp, UInt side_this, UInt side_other) except + nogil 

        # /// set an Id which allows unique identification of a compomer
        void setID(Size id) except + nogil  # wrap-doc:Sets an Id which allows unique identification of a compomer

        # /// return Id which allows unique identification of this compomer
        Size getID() except + nogil  # wrap-doc:Returns Id which allows unique identification of this compomer

        # /// left and right adducts of this compomer
        libcpp_vector[libcpp_map[String, Adduct] ] getComponent() except + nogil  # wrap-ignore
    
        # /// net charge of compomer (i.e. difference between left and right side of compomer)
        Int getNetCharge() except + nogil  # wrap-doc:Net charge of compomer (i.e. difference between left and right side of compomer)

        # /// mass of all contained adducts
        double getMass() except + nogil  # wrap-doc:Mass of all contained adducts

        # /// summed positive charges of contained adducts
        Int getPositiveCharges() except + nogil  # wrap-doc:Summed positive charges of contained adducts

        # /// summed negative charges of contained adducts
        Int getNegativeCharges() except + nogil  # wrap-doc:Summed negative charges of contained adducts

        # /// return log probability
        double getLogP() except + nogil  # wrap-doc:Returns the log probability

        # /// return log probability
        double getRTShift() except + nogil  # wrap-doc:Returns the log probability

        # /// get adducts with their abundance as compact string for both sides
        String getAdductsAsString() except + nogil  # wrap-doc:Get adducts with their abundance as compact string for both sides

        # /// get adducts with their abundance as compact string (amounts are absolute unless side=BOTH)
        # /// @param side Use LEFT for left, RIGHT for right
        String getAdductsAsString(UInt side) except + nogil  # wrap-doc:Get adducts with their abundance as compact string (amounts are absolute unless side=BOTH)

        # /// check if Compomer only contains a single adduct on side @p side
        bool isSingleAdduct(Adduct & a, UInt side) except + nogil  # wrap-doc:Check if Compomer only contains a single adduct on side @p side

        Compomer removeAdduct(Adduct & a) except + nogil  # wrap-doc:Remove ALL instances of the given adduct

        Compomer removeAdduct(Adduct & a, UInt side) except + nogil 

        StringList getLabels(UInt side) except + nogil  # wrap-doc:Returns the adduct labels from parameter(side) given. (LEFT or RIGHT)

        # void add(CompomerSide & add_side, UInt side) except + nogil 

cdef extern from "<OpenMS/DATASTRUCTURES/Compomer.h>" namespace "OpenMS::Compomer":
    # side of compomer (LEFT ^ substract; RIGHT ^ add)
    cdef enum SIDE: LEFT, RIGHT, BOTH
