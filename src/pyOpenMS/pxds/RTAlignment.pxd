from QCBase cimport QCBase
from FeatureMap cimport FeatureMap
from PeptideIdentificationList cimport PeptideIdentificationList
from TransformationDescription cimport TransformationDescription
from String cimport String


cdef extern from "<OpenMS/QC/RTAlignment.h>" namespace "OpenMS":

    cdef cppclass RTAlignment(QCBase):

        RTAlignment() except +

        void compute(FeatureMap& fm,
                     const TransformationDescription& trafo) const except +

        void compute(PeptideIdentificationList& ids,
                     const TransformationDescription& trafo) const except +

        const String& getName() const except +

        QCBase.Status requirements() const except +