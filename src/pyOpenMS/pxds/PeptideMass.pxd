from Types cimport *
from String cimport *
from QCBase cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/QC/PeptideMass.h>" namespace "OpenMS":
    
    cdef cppclass PeptideMass(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  QC metric calculating theoretical mass of a peptide sequence
        #  
        #  Each PeptideHit in the FeatureMap will be annotated with its 
        #  theoretical mass as metavalue 'mass'.

        PeptideMass() except + nogil 
        PeptideMass(PeptideMass &) except + nogil 

        void compute(FeatureMap & features) except + nogil 
            # wrap-doc:
            #  Sets the 'mass' metavalue to all PeptideHits by computing the theoretical mass
            #  
            #  :param features: FeatureMap with PeptideHits

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric
