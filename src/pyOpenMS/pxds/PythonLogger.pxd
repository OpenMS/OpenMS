from LogStream cimport *

# see ../extra_includes/python_logger.hpp for actual wrapped C++ code
cdef extern from "python_logger.hpp":

    cdef void redirectLoggingToPython() except +












