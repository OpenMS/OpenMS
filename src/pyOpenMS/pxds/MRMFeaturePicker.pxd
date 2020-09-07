from Param cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeaturePicker.h>" namespace "OpenMS":

    cdef cppclass MRMFeaturePicker:

        MRMFeaturePicker() nogil except +
        MRMFeaturePicker(MRMFeaturePicker &) nogil except +

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeaturePicker.h>" namespace "OpenMS::MRMFeaturePicker":

    cdef cppclass MRMFP_ComponentParams "OpenMS::MRMFeaturePicker::ComponentParams":

        MRMFP_ComponentParams() nogil except +
        MRMFP_ComponentParams(MRMFP_ComponentParams &) nogil except +

        String component_name
        String component_group_name
        Param params

    cdef cppclass MRMFP_ComponentGroupParams "OpenMS::MRMFeaturePicker::ComponentGroupParams":

        MRMFP_ComponentGroupParams() nogil except +
        MRMFP_ComponentGroupParams(MRMFP_ComponentGroupParams &) nogil except +

        String component_group_name
        Param params
