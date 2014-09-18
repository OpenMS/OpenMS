from Types cimport *
from smart_ptr cimport shared_ptr
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *
from LightTargetedExperiment cimport *
from AASequence cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>" namespace "OpenMS":

    cdef cppclass DIAScoring:
        DIAScoring() nogil except +
        # DIAScoring(DIAScoring) nogil except + #private

        void set_dia_parameters(double dia_extract_window, double dia_centroided,
                                double dia_byseries_intensity_min, double
                                dia_byseries_ppm_diff, double dia_nr_isotopes,
                                double dia_nr_charges) nogil except +

        bool dia_ms1_massdiff_score(double precursor_mz, OSSpectrumPtr spectrum,
                                    double& ppm_score) nogil except + 

        void dia_ms1_isotope_scores(double precursor_mz, OSSpectrumPtr spectrum, size_t charge_state, 
                                    double& isotope_corr, double& isotope_overlap) nogil except +

        # TODO automatically wrap 
        void dia_by_ion_score(OSSpectrumPtr spectrum, AASequence sequence, int charge, double & bseries_score, double & yseries_score) nogil except + # wrap-return:return(bseries_score,yseries_score) wrap-ignore

        # Dotproduct / Manhatten score with theoretical spectrum
        void score_with_isotopes(OSSpectrumPtr spectrum, libcpp_vector[LightTransition] transitions,
                                 double& dotprod, double& manhattan) nogil except +

