from Types cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from ILPDCWrapper cimport *
from DPosition cimport *
from MassExplainer cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.h>" namespace "OpenMS":
    
    cdef cppclass MetaboliteFeatureDeconvolution(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        MetaboliteFeatureDeconvolution() nogil except +
        MetaboliteFeatureDeconvolution(MetaboliteFeatureDeconvolution) nogil except +
        void compute(FeatureMap & fm_in, FeatureMap & fm_out, ConsensusMap & cons_map, ConsensusMap & cons_map_p) nogil except +

cdef extern from "<OpenMS/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.h>" namespace "OpenMS::MetaboliteFeatureDeconvolution":

    cdef enum CHARGEMODE_MFD "OpenMS::MetaboliteFeatureDeconvolution::CHARGEMODE":
        #wrap-attach:
        #    MetaboliteFeatureDeconvolution
        #wrap-instances:
        #    CHARGEMODE := CHARGEMODE_MFD

        QFROMFEATURE
        QHEURISTIC
        QALL

