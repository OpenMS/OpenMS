from Types cimport *
from String cimport *
from StringList cimport *

from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FORMAT/GNPSMGFFile.h>" namespace "OpenMS":

    cdef cppclass GNPSMGFFile(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler        

        GNPSMGFFile() except + nogil 
        GNPSMGFFile(GNPSMGFFile &) except + nogil 

        void store(const String& consensus_file_path, const StringList& mzml_file_paths, const String& out) except + nogil  # wrap-doc:Export consensus file from default workflow to GNPS MGF format
