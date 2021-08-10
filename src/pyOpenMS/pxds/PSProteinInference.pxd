from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.set cimport set as libcpp_set
from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from LPWrapper cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/PSProteinInference.h>" namespace "OpenMS":
    
    cdef cppclass PSProteinInference "OpenMS::PSProteinInference":
        PSProteinInference() nogil except + # wrap-doc:This class implements protein inference for the precursor ion selection strategies
        PSProteinInference(PSProteinInference &) nogil except + # compiler
        Size findMinimalProteinList(libcpp_vector[ PeptideIdentification ] & peptide_ids) nogil except +
        void calculateProteinProbabilities(libcpp_vector[ PeptideIdentification ] & ids) nogil except +
        double getProteinProbability(const String & acc) nogil except +
        bool isProteinInMinimalList(const String & acc) nogil except +
        Int getNumberOfProtIds(double protein_id_threshold) nogil except +
        # TODO nested STL
        # Int getNumberOfProtIdsPeptideRule(Int min_peptides, libcpp_map[ String, libcpp_set[ String ] ] & prot_id_counter) nogil except +
        SOLVER getSolver() nogil except +

