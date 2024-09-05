from InterpolationModel cimport *
from IsotopeDistribution cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/FEATUREFINDER/IsotopeModel.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeModel "OpenMS::IsotopeModel":
        # wrap-doc:
                #  Isotope distribution approximated using linear interpolation
                #  
                #  This models a smoothed (widened) distribution, i.e. can be used to sample actual raw peaks (depending on the points you query)
                #  If you only want the distribution (no widening), use either
                #  EmpiricalFormula::getIsotopeDistribution() // for a certain sum formula
                #  or
                #  IsotopeDistribution::estimateFromPeptideWeight (double average_weight)  // for averagine
                #  
                #  Peak widening is achieved by either a Gaussian or Lorentzian shape

        IsotopeModel() except + nogil 
        IsotopeModel(IsotopeModel &) except + nogil 
        UInt getCharge() except + nogil 
        void setOffset(double offset) except + nogil 
            # wrap-doc:
                #  Set the offset of the model
                #  
                #  The whole model will be shifted to the new offset without being computing all over
                #  This leaves a discrepancy which is minor in small shifts (i.e. shifting by one or two
                #  standard deviations) but can get significant otherwise. In that case use setParameters()
                #  which enforces a recomputation of the model

        double getOffset() except + nogil  # wrap-doc:Get the offset of the model
        EmpiricalFormula getFormula() except + nogil  # wrap-doc:Return the Averagine peptide formula (mass calculated from mean mass and charge -- use .setParameters() to set them)
        void setSamples(EmpiricalFormula &formula) except + nogil  # wrap-doc:Set sample/supporting points of interpolation
        double getCenter() except + nogil 
            # wrap-doc:
                #  Get the center of the Isotope model
                #  
                #  This is a m/z-value not necessarily the monoisotopic mass

        IsotopeDistribution  getIsotopeDistribution() except + nogil 
            # wrap-doc:
                #  Get the Isotope distribution (without widening) from the last setSamples() call
                #  
                #  Useful to determine the number of isotopes that the model contains and their position

        # BaseModel * create() except + nogil 
        
cdef extern from "<OpenMS/FEATUREFINDER/IsotopeModel.h>" namespace "OpenMS::IsotopeModel":
    cdef enum Averagines "OpenMS::IsotopeModel::Averagines":
        #wrap-attach:
        #   IsotopeModel
        C
        H
        N
        O
        S
        AVERAGINE_NUM

