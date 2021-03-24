from Types cimport *
from String cimport *
from smart_ptr cimport shared_ptr
from OpenSwathDataStructures cimport *
from LightTargetedExperiment cimport *
from AASequence cimport *
from DefaultParamHandler cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>" namespace "OpenMS":

    cdef cppclass DIAScoring(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        DIAScoring() nogil except +
        # DIAScoring(DIAScoring) nogil except + #private

        bool dia_ms1_massdiff_score(double precursor_mz, OSSpectrumPtr spectrum,
                                    double& ppm_score) nogil except + 

        void dia_ms1_isotope_scores_averagine(double precursor_mz, OSSpectrumPtr spectrum,
                                    double& isotope_corr, double& isotope_overlap, int charge_state) nogil except +

        void dia_ms1_isotope_scores(double precursor_mz, OSSpectrumPtr spectrum,
                                    double& isotope_corr, double& isotope_overlap, EmpiricalFormula& sum_formula) nogil except +
        # TODO automatically wrap 
        void dia_by_ion_score(OSSpectrumPtr spectrum, AASequence sequence, int charge, double & bseries_score, double & yseries_score) nogil except + # wrap-return:return(bseries_score,yseries_score) wrap-ignore

        # Dotproduct / Manhatten score with theoretical spectrum
        void score_with_isotopes(OSSpectrumPtr spectrum, libcpp_vector[LightTransition] transitions,
                                 double& dotprod, double& manhattan) nogil except +

