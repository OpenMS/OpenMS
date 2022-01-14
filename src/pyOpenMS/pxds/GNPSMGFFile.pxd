from Types cimport *
from String cimport *
from StringList cimport *

from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FORMAT/GNPSMGFFile.h>" namespace "OpenMS":

    cdef cppclass GNPSMGFFile(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler        

        GNPSMGFFile() nogil except +
        GNPSMGFFile(GNPSMGFFile &) nogil except +

        void run(String & consensus_file_path, StringList & mzml_file_paths, String & out) nogil except + # wrap-doc:Export consensus file from default workflow to GNPS MGF format
        
