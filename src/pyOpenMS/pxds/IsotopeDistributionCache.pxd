from Types cimport *
from FeatureFinderAlgorithmPickedHelperStructs cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/IsotopeDistributionCache.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeDistributionCache "OpenMS::IsotopeDistributionCache":
        IsotopeDistributionCache(IsotopeDistributionCache) nogil except + #wrap-ignore
        IsotopeDistributionCache(double max_mass, double mass_window_width, double intensity_percentage, double intensity_percentage_optional) nogil except +
        TheoreticalIsotopePattern  getIsotopeDistribution(double mass) nogil except +

