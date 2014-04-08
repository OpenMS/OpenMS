from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from FeatureMap cimport *
from Feature cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from PeptideIdentification cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from Types cimport *

cdef extern from "<OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>" namespace "OpenMS":

    cdef cppclass InternalCalibration(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        InternalCalibration()      nogil except +
        InternalCalibration(InternalCalibration)      nogil except + 
        void calibrateMapSpectrumwise(MSExperiment[Peak1D,ChromatogramPeak] & raw,
                                      MSExperiment[Peak1D,ChromatogramPeak] & calibrated,
                                      libcpp_vector[double] & ref_masses)      nogil except +
        void calibrateMapGlobally(MSExperiment[Peak1D,ChromatogramPeak] & raw,
                                  MSExperiment[Peak1D,ChromatogramPeak] & calibrated,
                                  libcpp_vector[double] & ref_masses,
                                  String trafo_filename)      nogil except +
        void calibrateMapGlobally(MSExperiment[Peak1D,ChromatogramPeak] & raw,
                                  MSExperiment[Peak1D,ChromatogramPeak] & calibrated,
                                  libcpp_vector[PeptideIdentification] & ref_ids,
                                  String trafo_filename)      nogil except +
        void calibrateMapGlobally(FeatureMap[Feature] & raw,
                                  FeatureMap[Feature] & calibrated,
                                  libcpp_vector[PeptideIdentification] & ref_ids,
                                  String trafo_filename)      nogil except +
        void calibrateMapGlobally(FeatureMap[Feature] & raw,
                                  FeatureMap[Feature] & calibrated,
                                  String trafo_filename)      nogil except +

