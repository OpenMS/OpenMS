from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair
# from libcpp cimport bool as libcpp_pair
from libcpp cimport bool 
from String cimport *
from DataValue cimport *
from CVTermList cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper":

    cdef cppclass RetentionTime:
        RetentionTime() nogil except +
        RetentionTime(RetentionTime) nogil except +
        String software_ref

    cdef cppclass Protein:

        Protein() nogil except +
        Protein(Protein) nogil except +

    cdef cppclass Peptide: # TODO add CVTermList

        Peptide() nogil except +
        Peptide(Peptide) nogil except +
        libcpp_vector[RetentionTime] rts
        String id
        libcpp_vector[String] protein_refs
        CVTermList evidence
        String sequence
        #std::vector<Modification> mods

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
