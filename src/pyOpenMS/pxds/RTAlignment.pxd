from QCBase cimport QCBase
from FeatureMap cimport FeatureMap
from PeptideIdentificationList cimport PeptideIdentificationList
from TransformationDescription cimport TransformationDescription
from String cimport String


cdef extern from "<OpenMS/QC/RTAlignment.h>" namespace "OpenMS":

    cdef cppclass RTAlignment(QCBase):

        RTAlignment() except +

        void compute(FeatureMap& fm,
                     const TransformationDescription& trafo) except + const 

        void compute(PeptideIdentificationList& ids,
                     const TransformationDescription& trafo)  except + const

        const String& getName() except + const

        QCBase.Status requirements() except + const
