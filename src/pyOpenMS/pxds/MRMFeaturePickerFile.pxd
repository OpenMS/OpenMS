from MRMFeaturePicker cimport *
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/FORMAT/MRMFeaturePickerFile.h>" namespace "OpenMS":

    cdef cppclass MRMFeaturePickerFile:
        # wrap-doc:
                #   _MRMFeaturePickerFile_ loads components and components groups parameters from a .csv file
                #   -----
                #   The structures defined in [MRMFeaturePicker](@ref MRMFeaturePicker) are used
                #   -----
                #   It is required that columns `component_name` and `component_group_name` are present.
                #   Lines whose `component_name`'s or `component_group_name`'s value is an empty string, will be skipped.
                #   The class supports the absence of information within other columns.
                #   -----
                #   A reduced example of the expected format (fewer columns are shown here):
                #   > component_name,component_group_name,TransitionGroupPicker:stop_after_feature,TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length
                #   > arg-L.arg-L_1.Heavy,arg-L,2,15
                #   > arg-L.arg-L_1.Light,arg-L,2,17
                #   > orn.orn_1.Heavy,orn,3,21
                #   > orn.orn_1.Light,orn,3,13

        MRMFeaturePickerFile() nogil except +
        MRMFeaturePickerFile(MRMFeaturePickerFile &) nogil except + # compiler

        void load(const String& filename, libcpp_vector[MRMFP_ComponentParams]& cp_list, libcpp_vector[MRMFP_ComponentGroupParams]& cgp_list) nogil except +
            # wrap-doc:
                #   Loads the file's data and saves it into vectors of `ComponentParams` and `ComponentGroupParams`
                #   -----
                #   The file is expected to contain at least two columns: `component_name` and `component_group_name`. Otherwise,
                #   an exception is thrown
                #   -----
                #   If a component group (identified by its name) is found multiple times, only the first one is saved
                #   -----
                #   :param filename: Path to the .csv input file
                #   :param cp_list: Component params are saved in this list
                #   :param cgp_list: Component Group params are saved in this list
                #   :raises:
                #     Exception: MissingInformation If the required columns are not found
                #   :raises:
                #     Exception: FileNotFound If input file is not found
