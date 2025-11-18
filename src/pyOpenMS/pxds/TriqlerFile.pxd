from String cimport *
from StringList cimport *
from ConsensusMap cimport *
from ExperimentalDesign cimport *

cdef extern from "<OpenMS/FORMAT/TriqlerFile.h>" namespace "OpenMS":

    cdef cppclass TriqlerFile:
        TriqlerFile() except + nogil
        TriqlerFile(TriqlerFile&) except + nogil
        void storeLFQ(const String& filename,
                      const ConsensusMap& consensus_map,
                      const ExperimentalDesign& design,
                      const StringList& reannotate_filenames,
                      const String& condition) except + nogil
