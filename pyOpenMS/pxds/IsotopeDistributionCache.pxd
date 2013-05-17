from Types cimport *
from FeatureFinderAlgorithmPickedHelperStructs cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/IsotopeDistributionCache.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeDistributionCache "OpenMS::IsotopeDistributionCache":
        IsotopeDistributionCache(IsotopeDistributionCache) nogil except + #wrap-ignore
        IsotopeDistributionCache(DoubleReal max_mass, DoubleReal mass_window_width, DoubleReal intensity_percentage, DoubleReal intensity_percentage_optional) nogil except +
        TheoreticalIsotopePattern  getIsotopeDistribution(DoubleReal mass) nogil except +

