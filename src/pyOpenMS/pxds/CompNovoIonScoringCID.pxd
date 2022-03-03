from PeptideIdentification cimport *
from DefaultParamHandler cimport *
from libcpp.map cimport map as cpp_map
# from ZhangSimilarityScore cimport *
# from MassDecomposition cimport *
# from MassDecompositionAlgorithm cimport *
# from CompNovoIonScoringBase cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringCID.h>" namespace "OpenMS":
    
    cdef cppclass CompNovoIonScoringCID "OpenMS::CompNovoIonScoringCID":
        CompNovoIonScoringCID() nogil except +
        CompNovoIonScoringCID(CompNovoIonScoringCID &) nogil except +
        # TODO OpenMS Map type
        # void scoreSpectrum(cpp_map[ double, IonScore ] & CID_ion_scores,
        #   MSSpectrum &CID_spec, double precursor_weight, Size charge) nogil except +

