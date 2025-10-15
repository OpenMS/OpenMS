from Types cimport *
from String cimport *
from QCBase cimport *
from FeatureMap cimport *
from TransformationDescription cimport *
from PeptideIdentificationList cimport *

cdef extern from "<OpenMS/QC/RTAlignment.h>" namespace "OpenMS":
    
    cdef cppclass RTAlignment(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  Take the original retention time before map alignment and use the alignment's 
        #  trafoXML for calculation of the new alignment retention times.
        #  
        #  Sets meta values "rt_raw" and "rt_align" in PeptideIdentifications of the 
        #  featureMap's PepIDs. It does not change the RT of the features.

        RTAlignment() except + nogil 
        RTAlignment(RTAlignment &) except + nogil 

        void compute(FeatureMap & fm, TransformationDescription & trafo) except + nogil 
            # wrap-doc:
            #  Calculates retention time after map alignment and sets meta values 
            #  "rt_raw" and "rt_align" in all PepIDs (on features and all unassigned PepIDs)
            #  
            #  :param fm: FeatureMap to receive the new metavalues
            #  :param trafo: Transformation information to get needed data from

        void compute(PeptideIdentificationList & ids, TransformationDescription & trafo) except + nogil 
            # wrap-doc:
            #  Calculates retention time after map alignment and sets meta values 
            #  "rt_raw" and "rt_align" in all PepIDs
            #  
            #  :param ids: PepIDs to receive the new metavalues
            #  :param trafo: Transformation information to get needed data from

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric
