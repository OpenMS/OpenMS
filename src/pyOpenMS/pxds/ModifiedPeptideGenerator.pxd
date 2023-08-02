from Types cimport *
from ProgressLogger cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from AASequence cimport *
from ResidueModification cimport *
from StringList cimport *

cdef extern from "<OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>" namespace "OpenMS":

    cdef cppclass ModifiedPeptideGenerator_MapToResidueType "OpenMS::ModifiedPeptideGenerator::MapToResidueType":

        ModifiedPeptideGenerator_MapToResidueType() except + nogil  # compiler
        ModifiedPeptideGenerator_MapToResidueType(ModifiedPeptideGenerator_MapToResidueType &) except + nogil  # compiler

    cdef cppclass ModifiedPeptideGenerator:
        # wrap-doc:
        #  Generates modified peptides/proteins.

        ModifiedPeptideGenerator() except + nogil  # compiler
        ModifiedPeptideGenerator(ModifiedPeptideGenerator &) except + nogil  # compiler

        @staticmethod
        ModifiedPeptideGenerator_MapToResidueType getModifications(const StringList& modNames) except + nogil 

        @staticmethod
        void applyFixedModifications(const ModifiedPeptideGenerator_MapToResidueType& fixed_mods, AASequence& peptide) except + nogil 

        @staticmethod
        void applyVariableModifications(const ModifiedPeptideGenerator_MapToResidueType& var_mods, const AASequence& peptide, Size max_variable_mods_per_peptide, libcpp_vector[AASequence]& all_modified_peptides, bool keep_original) except + nogil 
