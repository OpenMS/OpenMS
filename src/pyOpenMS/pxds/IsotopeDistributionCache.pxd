from Types cimport *
from FeatureFinderAlgorithmPickedHelperStructs cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/IsotopeDistributionCache.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeDistributionCache "OpenMS::IsotopeDistributionCache":
        IsotopeDistributionCache(double max_mass, double mass_window_width, double intensity_percentage, double intensity_percentage_optional) except + nogil 
        IsotopeDistributionCache(IsotopeDistributionCache &) except + nogil 
        TheoreticalIsotopePattern  getIsotopeDistribution(double mass) except + nogil  # wrap-doc:Returns the isotope distribution for a certain mass window
