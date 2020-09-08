from Types cimport *
from ProgressLogger cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from AASequence cimport *
from ResidueModification cimport *
from StringList cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>" namespace "OpenMS":

    cdef cppclass ModifiedPeptideGenerator:

        ModifiedPeptideGenerator() nogil except +
        ModifiedPeptideGenerator(ModifiedPeptideGenerator) nogil except + 

    cdef cppclass ModifiedPeptideGenerator_MapToResidueType "OpenMS::ModifiedPeptideGenerator::MapToResidueType":

        ModifiedPeptideGenerator_MapToResidueType() nogil except +
        ModifiedPeptideGenerator_MapToResidueType(ModifiedPeptideGenerator_MapToResidueType) nogil except + #wrap-ignore


    ## wrap static methods
    cdef extern from "<OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>" namespace "OpenMS::ModifiedPeptideGenerator":
        ModifiedPeptideGenerator_MapToResidueType getModifications(const StringList& modNames) nogil except +

        void applyFixedModifications(const ModifiedPeptideGenerator_MapToResidueType& fixed_mods, AASequence& peptide) nogil except +

        void applyVariableModifications(const ModifiedPeptideGenerator_MapToResidueType& var_mods, 
          const AASequence& peptide, 
          Size max_variable_mods_per_peptide, 
          libcpp_vector[AASequence]& all_modified_peptides, 
          bool keep_original) nogil except +

