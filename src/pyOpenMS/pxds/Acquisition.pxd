from libcpp cimport bool
from MetaInfoInterface cimport *
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/METADATA/Acquisition.h>" namespace "OpenMS":
	
	cdef cppclass Acquisition(MetaInfoInterface):
		# wrap-inherits:
		#  MetaInfoInterface

		Acquisition() except + nogil 
		Acquisition(Acquisition &) except + nogil 

		bool operator==(Acquisition &rhs) except + nogil 
		bool operator!=(Acquisition &rhs) except + nogil 
		String getIdentifier() except + nogil 
		void setIdentifier(const String &identifier) except + nogil 
