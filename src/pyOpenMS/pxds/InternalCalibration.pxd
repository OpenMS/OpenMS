from libcpp.vector cimport vector as libcpp_vector
from FeatureMap cimport *
from MZTrafoModel cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from PeptideIdentification cimport *
from ProgressLogger cimport *
from Types cimport *

cdef extern from "<OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>" namespace "OpenMS":

    cdef cppclass InternalCalibration(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        InternalCalibration()      nogil except +
        InternalCalibration(InternalCalibration) nogil except + 
        
        Size fillCalibrants(MSExperiment[Peak1D,ChromatogramPeak], libcpp_vector[InternalCalibration::LockMass], double tol_ppm, bool, bool, bool) nogil except +
        Size fillCalibrants(FeatureMap, double) nogil except +
        Size fillCalibrants(libcpp_vector[PeptideIdentification], double) nogil except +
        CalibrationData getCalibrationPoints() nogil except +
        bool calibrate(MSExperiment[Peak1D,ChromatogramPeak], libcpp_vector[int], MZTrafoModel::MODELTYPE, double, bool, double, double, String, String) nogil except +
                                     
                      
cdef extern from "<OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>" namespace "OpenMS::InternalCalibration":
    
    # static members               
    static void applyTransformation(MSExperiment<>::SpectrumType& spec, const MZTrafoModel& trafo) nogil except +
    static void applyTransformation(MSExperiment<>& exp, const IntList& target_mslvl, const MZTrafoModel& trafo) nogil except +