from Param cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeaturePicker.h>" namespace "OpenMS":

    cdef cppclass MRMFeaturePicker:

        MRMFeaturePicker() except + nogil 
        MRMFeaturePicker(MRMFeaturePicker &) except + nogil  # compiler

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeaturePicker.h>" namespace "OpenMS::MRMFeaturePicker":

    cdef cppclass MRMFP_ComponentParams "OpenMS::MRMFeaturePicker::ComponentParams":

        MRMFP_ComponentParams() except + nogil 
        MRMFP_ComponentParams(MRMFP_ComponentParams &) except + nogil 

        String component_name
        String component_group_name
        Param params

    cdef cppclass MRMFP_ComponentGroupParams "OpenMS::MRMFeaturePicker::ComponentGroupParams":

        MRMFP_ComponentGroupParams() except + nogil 
        MRMFP_ComponentGroupParams(MRMFP_ComponentGroupParams &) except + nogil 

        String component_group_name
        Param params
