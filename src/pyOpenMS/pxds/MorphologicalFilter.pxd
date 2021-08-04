from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>" namespace "OpenMS":

    cdef cppclass MorphologicalFilter(DefaultParamHandler,ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        MorphologicalFilter() nogil except +
        # private
        MorphologicalFilter(MorphologicalFilter) nogil except + # wrap-ignore

        void filter(MSSpectrum & spectrum)      nogil except +
            # wrap-doc:
                #   Applies the morphological filtering operation to an MSSpectrum
                #   -----
                #   If the size of the structuring element is given in 'Thomson', the number of data points for
                #   the structuring element is computed as follows:
                #
                #       - The data points are assumed to be uniformly spaced.  We compute the
                #           average spacing from the position of the first and the last peak and the
                #           total number of peaks in the input range
                #       - The number of data points in the structuring element is computed
                #           from struc_size and the average spacing, and rounded up to an odd
                #           number
                
        void filterExperiment(MSExperiment & exp)      nogil except +
            # wrap-doc:
                #   Applies the morphological filtering operation to an MSExperiment
                #   -----
                #   The size of the structuring element is computed for each spectrum individually, if it is given in 'Thomson'
                #   See the filtering method for MSSpectrum for details
