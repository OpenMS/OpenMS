from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.set cimport set as libcpp_set
from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from LPWrapper cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/PSProteinInference.h>" namespace "OpenMS":
    
    cdef cppclass PSProteinInference "OpenMS::PSProteinInference":
        # wrap-doc:
        # This class implements protein inference for the precursor ion selection strategies
        PSProteinInference() except + nogil 
        PSProteinInference(PSProteinInference &) except + nogil  # compiler
        Size findMinimalProteinList(libcpp_vector[ PeptideIdentification ] & peptide_ids) except + nogil 
        void calculateProteinProbabilities(libcpp_vector[ PeptideIdentification ] & ids) except + nogil 
        double getProteinProbability(const String & acc) except + nogil 
        bool isProteinInMinimalList(const String & acc) except + nogil 
        Int getNumberOfProtIds(double protein_id_threshold) except + nogil 
        # TODO nested STL
        # Int getNumberOfProtIdsPeptideRule(Int min_peptides, libcpp_map[ String, libcpp_set[ String ] ] & prot_id_counter) except + nogil 
        SOLVER getSolver() except + nogil 

