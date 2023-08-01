from PeptideIdentification cimport *
from DefaultParamHandler cimport *
from Types cimport *
# from ZhangSimilarityScore cimport *
# from MassDecomposition cimport *
# from MassDecompositionAlgorithm cimport *
# from CompNovoIonScoringBase cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringCID.h>" namespace "OpenMS":
    
    cdef cppclass CompNovoIonScoringCID "OpenMS::CompNovoIonScoringCID":
        CompNovoIonScoringCID() except + nogil 
        CompNovoIonScoringCID(CompNovoIonScoringCID &) except + nogil 
        # TODO OpenMS Map type
        # void scoreSpectrum(libcpp_map[ double, IonScore ] & CID_ion_scores,
        #  MSSpectrum &CID_spec, double precursor_weight, Size charge) except + nogil 

