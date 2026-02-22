from QCBase cimport QCBase
from FeatureMap cimport FeatureMap
from String cimport String


cdef extern from "<OpenMS/QC/FWHM.h>" namespace "OpenMS":

    cdef cppclass FWHM(QCBase):

        FWHM() except +

        void compute(FeatureMap& features) except +

        const String& getName() const except +

        QCBase.Status requirements() const except +