from PeptideIdentification cimport *
from DefaultParamHandler cimport *
from Map cimport *
# from ZhangSimilarityScore cimport *
# from MassDecomposition cimport *
# from MassDecompositionAlgorithm cimport *
# from CompNovoIonScoringBase cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringCID.h>" namespace "OpenMS":
    
    cdef cppclass CompNovoIonScoringCID "OpenMS::CompNovoIonScoringCID":
        CompNovoIonScoringCID() nogil except +
        CompNovoIonScoringCID(CompNovoIonScoringCID) nogil except +
        # void scoreSpectrum(Map[ DoubleReal, IonScore ] & CID_ion_scores, 
        #   MSSpectrum[Peak1D] &CID_spec, DoubleReal precursor_weight, Size charge) nogil except +

