from libcpp cimport bool
from MetaInfoInterface cimport *
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/Acquisition.h>" namespace "OpenMS":
	
	cdef cppclass Acquisition(MetaInfoInterface):
		# wrap-inherits:
		#    MetaInfoInterface

		Acquisition() nogil except +
		Acquisition(Acquisition) nogil except +

		bool operator==(Acquisition &rhs) nogil except +
		bool operator!=(Acquisition &rhs) nogil except +
		String getIdentifier() nogil except +
		void setIdentifier(const String &identifier) nogil except +

		void getKeys(libcpp_vector[String] & keys)
		DataValue getMetaValue(String) nogil except +
		void setMetaValue(String, DataValue) nogil except +
		bool metaValueExists(String) nogil except +
		void removeMetaValue(String) nogil except +	