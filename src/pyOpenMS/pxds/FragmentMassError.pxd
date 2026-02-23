from libcpp.vector cimport vector

from QCBase cimport QCBase
from FeatureMap cimport FeatureMap
from MSExperiment cimport MSExperiment
from PeptideIdentificationList cimport PeptideIdentificationList
from ProteinIdentification cimport ProteinIdentification
from Types cimport UInt32
from String cimport String


cdef extern from "<OpenMS/QC/FragmentMassError.h>" namespace "OpenMS":

    cdef cppclass FragmentMassError(QCBase):
        # Nested struct Statistics
        cppclass Statistics:
            double average_ppm
            double variance_ppm

        
        # Enum ToleranceUnit
        
        cdef enum ToleranceUnit:
            AUTO
            PPM
            DA

        FragmentMassError() except +

        
        # compute overload 1
        
        void compute(FeatureMap& fmap,
                     const MSExperiment& exp,
                     const QCBase.SpectraMap& map_to_spectrum,
                     ToleranceUnit tolerance_unit = AUTO,
                     double tolerance = 20) except +

       
        # compute overload 2
        
        void compute(PeptideIdentificationList& pep_ids,
                     const ProteinIdentification.SearchParameters& search_params,
                     const MSExperiment& exp,
                     const QCBase.SpectraMap& map_to_spectrum,
                     ToleranceUnit tolerance_unit = AUTO,
                     double tolerance = 20) except +

        const String& getName() except + const

        const vector[Statistics]& getResults() except + const

        QCBase.Status requirements() except + const