from Types cimport *
from String cimport *
from smart_ptr cimport shared_ptr
from OpenSwathDataStructures cimport *
from LightTargetedExperiment cimport *
from AASequence cimport *
from DefaultParamHandler cimport *
from EmpiricalFormula cimport *
from libcpp.vector cimport vector as libcpp_vector
from RangeManager cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>" namespace "OpenMS":

    cdef cppclass DIAScoring(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        DIAScoring() except + nogil 
        # private
        DIAScoring(DIAScoring) except + nogil  # wrap-ignore

        bool dia_ms1_massdiff_score(double precursor_mz, libcpp_vector[shared_ptr[OSSpectrum] ] spectrum,
                                    RangeMobility& im_range, double& ppm_score) except + nogil

        void dia_ms1_isotope_scores_averagine(double precursor_mz, libcpp_vector[shared_ptr[OSSpectrum] ] spectrum, int charge_state, RangeMobility& im_range,
                double& isotope_corr, double& isotope_overlap) except + nogil

        void dia_ms1_isotope_scores(double precursor_mz, libcpp_vector[shared_ptr[OSSpectrum] ] spectrum, RangeMobility& im_range,
                double& isotope_corr, double& isotope_overlap, EmpiricalFormula& sum_formula) except + nogil
        # TODO automatically wrap 
        void dia_by_ion_score( libcpp_vector[shared_ptr[OSSpectrum] ] spectrum, AASequence sequence, int charge, RangeMobility& im_range, double & bseries_score, double & yseries_score) except + nogil # wrap-return:return(bseries_score,yseries_score) wrap-ignore

        # Dotproduct / Manhatten score with theoretical spectrum
        void score_with_isotopes( libcpp_vector[shared_ptr[OSSpectrum] ] spectrum, libcpp_vector[LightTransition] transitions, RangeMobility& im_range,
                                 double& dotprod, double& manhattan) except + nogil
