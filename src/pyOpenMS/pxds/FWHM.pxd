from Types cimport *
from String cimport *
from QCBase cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/QC/FWHM.h>" namespace "OpenMS":
    
    cdef cppclass FWHM(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  QC metric calculating Full Width at Half Maximum (FWHM)
        #  
        #  The metric moves FWHM metavalues from the feature to all its 
        #  PeptideIdentifications (since that's where mzTab takes it from 
        #  if we want to preserve raw file origin).

        FWHM() except + nogil 
        FWHM(FWHM &) except + nogil 

        void compute(FeatureMap & features) except + nogil 
            # wrap-doc:
            #  Moves FWHM metavalues from the feature to all its PeptideIdentifications
            #  
            #  A warning is issued on the commandline if a feature does not have 
            #  either 'FWHM' or 'model_FWHM' as metavalue.
            #  
            #  :param features: FeatureMap with features which have metavalue 'FWHM' or 'model_FWHM'

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric
