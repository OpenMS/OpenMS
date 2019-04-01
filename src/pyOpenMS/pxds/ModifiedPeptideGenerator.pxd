from ProgressLogger cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from AASequence cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>" namespace "OpenMS":

    cdef cppclass ModifiedPeptideGenerator:

        ModifiedPeptideGenerator() nogil except +
        ModifiedPeptideGenerator(ModifiedPeptideGenerator) nogil except + 

    ## wrap static methods
    cdef extern from "<OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>" namespace "OpenMS::ModifiedPeptideGenerator":
        void applyFixedModifications(const libcpp_vector[const ResidueModification*]& fixed_mods, AASequence& peptide);

        void applyVariableModifications(const libcpp_vector[const ResidueModification*]& var_mods, 
          const AASequence& peptide, 
          Size max_variable_mods_per_peptide, 
          libcpp_vector[AASequence]& all_modified_peptides, 
          bool keep_original) nogil except +

