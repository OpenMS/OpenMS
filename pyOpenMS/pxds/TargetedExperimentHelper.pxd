from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from String cimport *
from CVTerm cimport *
from Map    cimport *
from DataValue cimport *
from CVTermList cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper":

    cdef cppclass RetentionTime:

        RetentionTime() nogil except +
        RetentionTime(RetentionTime) nogil except +
        String software_ref

        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +

        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms,
                            String accession
                           ) nogil except +

        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map
                           ) nogil except +

        Map[String, libcpp_vector[CVTerm] ] getCVTerms()
        void addCVTerm(CVTerm & term)                   nogil except +

        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

    cdef cppclass Protein:

        Protein() nogil except +
        Protein(Protein) nogil except +

        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +

        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms,
                            String accession
                           ) nogil except +

        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map
                           ) nogil except +

        Map[String, libcpp_vector[CVTerm] ] getCVTerms()
        void addCVTerm(CVTerm & term)                   nogil except +

        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

    cdef cppclass Peptide:

        Peptide() nogil except +
        Peptide(Peptide) nogil except +

        # members
        libcpp_vector[RetentionTime] rts
        String id
        libcpp_vector[String] protein_refs
        CVTermList evidence
        String sequence
        libcpp_vector[Modification] mods

        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +

        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms,
                            String accession
                           ) nogil except +

        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map
                           ) nogil except +

        Map[String, libcpp_vector[CVTerm] ] getCVTerms()
        void addCVTerm(CVTerm & term)                   nogil except +

        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

        void setChargeState(int charge) nogil except +
        int getChargeState() nogil except +
        void setPeptideGroupLabel(String label) nogil except +
        String getPeptideGroupLabel() nogil except +
        double getRetentionTime() nogil except +

        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys)
        void getKeys(libcpp_vector[unsigned int] & keys)
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +

# no support for nested classes yet in Cython
cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper::Peptide":

    cdef cppclass Modification:

        Modification() nogil except +
        Modification(Modification) nogil except +

        # members
        DoubleReal avg_mass_delta
        int location
        DoubleReal mono_mass_delta

