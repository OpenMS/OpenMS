from ProgressLogger cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>" namespace "OpenMS":

    cdef cppclass ModifiedPeptideGenerator:

        ModifiedPeptideGenerator() nogil except +
        ModifiedPeptideGenerator(ModifiedPeptideGenerator) nogil except + 

        # TODO : iterator as arg
        # applyFixedModifications

        # TODO : iterator as arg
        # applyVariableModifications
