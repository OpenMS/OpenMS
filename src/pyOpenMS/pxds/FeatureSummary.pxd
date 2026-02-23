from libcpp.vector cimport vector

from QCBase cimport QCBase
from FeatureMap cimport FeatureMap
from Types cimport UInt
from String cimport String


cdef extern from "<OpenMS/QC/FeatureSummary.h>" namespace "OpenMS":

    cdef cppclass FeatureSummary(QCBase):

        cppclass Result:
            UInt feature_count
            float rt_shift_mean

        FeatureSummary() except +

        Result compute(const FeatureMap& feature_map) except +

        const String& getName() except + const

        QCBase.Status requirements() except + const