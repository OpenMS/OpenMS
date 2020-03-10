from MRMFeaturePicker cimport *
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/FORMAT/MRMFeaturePickerFile.h>" namespace "OpenMS":

    cdef cppclass MRMFeaturePickerFile:

        MRMFeaturePickerFile() nogil except +
        MRMFeaturePickerFile(MRMFeaturePickerFile &) nogil except +

        void load(const String& filename, libcpp_vector[MRMFP_ComponentParams]& cp_list, libcpp_vector[MRMFP_ComponentGroupParams]& cgp_list) nogil except +
