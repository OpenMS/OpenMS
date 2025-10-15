from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from Types cimport *
from String cimport *
from QCBase cimport *
from FeatureMap cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/QC/MissedCleavages.h>" namespace "OpenMS":
    
    cdef cppclass MissedCleavages(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  This class is a metric for the QualityControl TOPP Tool.
        #  
        #  This class counts the number of missed cleavages per PeptideIdentification 
        #  given a FeatureMap and returns an agglomeration statistic (observed counts).
        #  Additionally the PeptideHits in the FeatureMap are augmented with meta information.

        MissedCleavages() except + nogil 
        MissedCleavages(MissedCleavages &) except + nogil 

        void compute(FeatureMap & fmap, SearchParameters & digestor, UInt32 max_mc) except + nogil 
            # wrap-doc:
            #  Counts the number of missed cleavages per PeptideIdentification.
            #  
            #  The result is a key/value map: #missed_cleavages --> counts
            #  
            #  :param fmap: Input FeatureMap
            #  :param digestor: Search parameters with digestion enzyme
            #  :param max_mc: Maximum number of missed cleavages to count

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric

        libcpp_vector[libcpp_map[UInt32, UInt32]] getResults() except + nogil  # wrap-doc:Returns results
