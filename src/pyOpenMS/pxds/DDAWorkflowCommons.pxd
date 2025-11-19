from libcpp.map cimport map as libcpp_map
from libcpp.string cimport string as libcpp_string
from String cimport *
from StringList cimport *
from MSExperiment cimport *
from FeatureMap cimport *
from PeptideIdentificationList cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/DDAWorkflowCommons.h>" namespace "OpenMS":

    cdef cppclass DDAWorkflowCommons:
        DDAWorkflowCommons() except + nogil
        DDAWorkflowCommons(DDAWorkflowCommons&) except + nogil

        libcpp_map[libcpp_string, libcpp_string] mapMzML2Ids(StringList& in_files, StringList& in_id_files) except + nogil
        libcpp_map[libcpp_string, libcpp_string] mapId2MzMLs(const libcpp_map[libcpp_string, libcpp_string]& mz_to_id) except + nogil
        double estimateMedianChromatographicFWHM(MSExperiment& ms_centroided) except + nogil
        void recalibrateMS1(MSExperiment& ms_centroided,
                           PeptideIdentificationList& peptide_ids,
                           const String& id_file_abs_path) except + nogil
        void calculateSeeds(const MSExperiment& ms_centroided,
                            double intensity_threshold,
                            FeatureMap& seeds,
                            double median_fwhm,
                            Size charge_min,
                            Size charge_max) except + nogil
