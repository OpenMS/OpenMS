from PeptideIdentification cimport *
from DefaultParamHandler cimport *
from libcpp.map cimport map as libcpp_map
# from ZhangSimilarityScore cimport *
# from MassDecomposition cimport *
# from MassDecompositionAlgorithm cimport *
# from CompNovoIonScoringBase cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/CompNovoIonScoring.h>" namespace "OpenMS":
    
    cdef cppclass CompNovoIonScoring "OpenMS::CompNovoIonScoring":
        CompNovoIonScoring() nogil except +
        CompNovoIonScoring(CompNovoIonScoring &) nogil except +
        # TODO -> replace type ... 
        # TODO OpenMS Map type
        # void scoreSpectra(libcpp_map[ double, IonScore ] &CID_ion_scores, MSSpectrum & CID_spec,
        #  MSSpectrum &ETD_spec, double precursor_weight, Size charge) nogil except +

