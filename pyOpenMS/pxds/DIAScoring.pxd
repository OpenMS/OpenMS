from smart_ptr cimport shared_ptr
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *
from AASequence cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>" namespace "OpenMS":

    cdef cppclass DIAScoring:
        DIAScoring() nogil except +
        # DIAScoring(DIAScoring) nogil except + #private

        void set_dia_parameters(double dia_extract_window, double dia_centroided,
                                double dia_byseries_intensity_min, double dia_byseries_ppm_diff, double dia_nr_isotopes, double dia_nr_charges) nogil except +

        # TODO automatically wrap 
        void dia_by_ion_score(shared_ptr[Spectrum] spectrum, AASequence sequence, int charge, double & bseries_score, double & yseries_score) nogil except + # wrap-return:return(bseries_score,yseries_score) wrap-ignore
