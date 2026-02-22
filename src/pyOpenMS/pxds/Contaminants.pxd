from libcpp.vector cimport vector
from libcpp.utility cimport pair

from QCBase cimport QCBase
from FeatureMap cimport FeatureMap
from FASTAFile cimport FASTAFile
from String cimport String
from Types cimport Int64
from PeptideHit cimport PeptideHit


cdef extern from "<OpenMS/QC/Contaminants.h>" namespace "OpenMS":

    cdef cppclass Contaminants(QCBase):

        # nested struct
        cdef cppclass ContaminantsSummary:
            double assigned_contaminants_ratio
            double unassigned_contaminants_ratio
            double all_contaminants_ratio
            double assigned_contaminants_intensity_ratio
            pair[Int64, Int64] empty_features

        Contaminants() except +

        void compute(FeatureMap& features,
                     const vector[FASTAFile.FASTAEntry]& contaminants) except +

        const String& getName() const except +

        const vector[ContaminantsSummary]& getResults() except +

        QCBase.Status requirements() const except +