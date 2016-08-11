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
        

        Size fillCalibrants(MSExperiment[Peak1D,ChromatogramPeak],
                            libcpp_vector[InternalCalibration_LockMass], 
                            double tol_ppm,
                            bool lock_require_mono,
                            bool lock_require_iso,
                            CalibrationData& failed_lock_masses,
                            bool verbose) nogil except +
        Size fillCalibrants(FeatureMap, double) nogil except +
        Size fillCalibrants(libcpp_vector[PeptideIdentification], double) nogil except +
        CalibrationData getCalibrationPoints() nogil except +
        bool calibrate(MSExperiment[Peak1D,ChromatogramPeak],
                       libcpp_vector[int], MZTrafoModel_MODELTYPE, double,
                       bool, double, double, String, String) nogil except +
                                     
## wrap static methods
cdef extern from "<OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>" namespace "OpenMS::InternalCalibration":
    
    void applyTransformation(MSSpectrum[Peak1D] & spec, 
                                    MZTrafoModel & trafo) nogil except + # wrap-attach:InternalCalibration
    void applyTransformation(MSExperiment[Peak1D, ChromatogramPeak] & exp, 
                                    IntList& target_mslvl, MZTrafoModel& trafo) nogil except + # wrap-attach:InternalCalibration

cdef extern from "<OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>" namespace "OpenMS::InternalCalibration":
    
    cdef cppclass InternalCalibration_LockMass "OpenMS::InternalCalibration::LockMass":
        #InternalCalibration_LockMass() nogil except +
        InternalCalibration_LockMass(InternalCalibration_LockMass) nogil except + # wrap-ignore
        InternalCalibration_LockMass(double mz_, int lvl_, int charge_) nogil except + 
        double mz # m/z of the lock mass (incl. adducts)
        unsigned int ms_level # MS level where it occurs
        int charge # charge of the ion (to find isotopes)

