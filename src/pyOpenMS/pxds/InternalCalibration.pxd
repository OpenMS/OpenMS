from libcpp.vector cimport vector as libcpp_vector
from FeatureMap cimport *
from MZTrafoModel cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from PeptideIdentification cimport *
from ProgressLogger cimport *
from Types cimport *
from Precursor cimport *

cdef extern from "<OpenMS/PROCESSING/CALIBRATION/InternalCalibration.h>" namespace "OpenMS":

    cdef cppclass InternalCalibration(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger

        InternalCalibration()      except + nogil  # wrap-doc:A mass recalibration method using linear/quadratic interpolation (robust/weighted) of given reference masses
        InternalCalibration(InternalCalibration &) except + nogil  # compiler


        Size fillCalibrants(MSExperiment,
                            libcpp_vector[InternalCalibration_LockMass],
                            double tol_ppm,
                            bool lock_require_mono,
                            bool lock_require_iso,
                            CalibrationData& failed_lock_masses,
                            bool verbose) except + nogil 
            # wrap-doc:
                #  Extract calibrants from Raw data (mzML)\n
                #  
                #  Lock masses are searched in each spectrum and added to the internal calibrant database\n
                #  
                #  Filters can be used to exclude spurious peaks, i.e. require the calibrant peak to be monoisotopic or
                #  to have a +1 isotope (should not be used for very low abundant calibrants)
                #  If a calibrant is not found, it is added to a 'failed_lock_masses' database which is returned and not stored internally.
                #  The intensity of the peaks describe the reason for failed detection: 0.0 - peak not found with the given ppm tolerance;
                #  1.0 - peak is not monoisotopic (can only occur if 'lock_require_mono' is true)
                #  2.0 - peak has no +1 isotope (can only occur if 'lock_require_iso' is true)
                #  
                #  
                #  :param exp: Peak map containing the lock masses
                #  :param ref_masses: List of lock masses
                #  :param tol_ppm: Search window for lock masses in 'exp'
                #  :param lock_require_mono: Require that a lock mass is the monoisotopic peak (i.e. not an isotope peak) -- lock mass is rejected otherwise
                #  :param lock_require_iso: Require that a lock mass has isotope peaks to its right -- lock mass is rejected otherwise
                #  :param failed_lock_masses: Set of calibration masses which were not found, i.e. their expected m/z and RT positions
                #  :param verbose: Print information on 'lock_require_XXX' matches during search
                #  :return: Number of calibration masses found

        Size fillCalibrants(FeatureMap, double) except + nogil 
            # wrap-doc:
                #  Extract calibrants from identifications\n
                #  
                #  Extracts only the first hit from the first peptide identification of each feature
                #  Hits are sorted beforehand
                #  Ambiguities should be resolved before, e.g. using IDFilter
                #  RT and m/z are taken from the features, not from the identifications (for an exception see below)!\n
                #  
                #  Unassigned peptide identifications are also taken into account!
                #  RT and m/z are naturally taken from the IDs, since to feature is assigned
                #  If you do not want these IDs, remove them from the feature map before calling this function\n
                #  
                #  A filtering step is done in the m/z dimension using 'tol_ppm'
                #  Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
                #  larger outliers are removed before accepting an ID as calibrant
                #  
                #  
                #  :param fm: FeatureMap with peptide identifications
                #  :param tol_ppm: Only accept ID's whose theoretical mass deviates at most this much from annotated
                #  :return: Number of calibration masses found

        Size fillCalibrants(libcpp_vector[PeptideIdentification], double) except + nogil 
            # wrap-doc:
                #  Extract calibrants from identifications\n
                #  
                #  Extracts only the first hit from each peptide identification
                #  Hits are sorted beforehand
                #  Ambiguities should be resolved before, e.g. using IDFilter\n
                #  
                #  Unassigned peptide identifications are also taken into account!
                #  RT and m/z are naturally taken from the IDs, since to feature is assigned
                #  If you do not want these IDs, remove them from the feature map before calling this function\n
                #  
                #  A filtering step is done in the m/z dimension using 'tol_ppm'
                #  Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
                #  larger outliers are removed before accepting an ID as calibrant
                #  
                #  
                #  :param pep_ids: Peptide ids (e.g. from an idXML file)
                #  :param tol_ppm: Only accept ID's whose theoretical mass deviates at most this much from annotated
                #  :return: Number of calibration masses found

        CalibrationData getCalibrationPoints() except + nogil 
            # wrap-doc:
                #  Get container of calibration points\n
                #  
                #  Filled using fillCalibrants() methods
                #  
                #  
                #  :return: Container of calibration points

        bool calibrate(MSExperiment,
                       libcpp_vector[int],
                       MZTrafoModel_MODELTYPE,
                       double rt_chunk,
                       bool use_RANSAC,
                       double post_ppm_median,
                       double post_ppm_MAD,
                       String file_models,
                       String file_models_plot,
                       String file_residuals,
                       String file_residuals_plot,
                       String rscript_executable) except + nogil 
            # wrap-doc:
                #  Apply calibration to data\n
                #  
                #  For each spectrum, a calibration model will be computed and applied.
                #  Make sure to call fillCalibrants() before, so a model can be created.\n
                #  
                #  The MSExperiment will be sorted by RT and m/z if unsorted.
                #  
                #  
                #  :param exp: MSExperiment holding the Raw data to calibrate
                #  :param target_mslvl: MS-levels where calibration should be applied to
                #  :param model_type: Linear or quadratic model; select based on your instrument
                #  :param rt_chunk: RT-window size (one-sided) of calibration points to collect around each spectrum. Set to negative values, to build one global model instead.
                #  :param use_RANSAC: Remove outliers before fitting a model?!
                #  :param post_ppm_median: The median ppm error of the calibrants must be at least this good after calibration; otherwise this method returns false(fail)
                #  :param post_ppm_MAD: The median absolute deviation of the calibrants must be at least this good after calibration; otherwise this method returns false(fail)
                #  :param file_models: Output CSV filename, where model parameters are written to (pass empty string to skip)
                #  :param file_models_plot: Output PNG image model parameters (pass empty string to skip)
                #  :param file_residuals: Output CSV filename, where ppm errors of calibrants before and after model fitting parameters are written to (pass empty string to skip)
                #  :param file_residuals_plot: Output PNG image of the ppm errors of calibrants (pass empty string to skip)
                #  :param rscript_executable: Full path to the Rscript executable
                #  :return: true upon successful calibration


## wrap static methods
cdef extern from "<OpenMS/PROCESSING/CALIBRATION/InternalCalibration.h>" namespace "OpenMS::InternalCalibration":

    void applyTransformation(libcpp_vector[Precursor]& pcs,
                             MZTrafoModel& trafo) except + nogil  # wrap-attach:InternalCalibration
    void applyTransformation(MSSpectrum & spec, IntList& target_mslvl,
                             MZTrafoModel & trafo) except + nogil  # wrap-attach:InternalCalibration
    void applyTransformation(MSExperiment & exp,
                             IntList& target_mslvl, MZTrafoModel& trafo) except + nogil  # wrap-attach:InternalCalibration

cdef extern from "<OpenMS/PROCESSING/CALIBRATION/InternalCalibration.h>" namespace "OpenMS::InternalCalibration":

    cdef cppclass InternalCalibration_LockMass "OpenMS::InternalCalibration::LockMass":
        InternalCalibration_LockMass(double mz_, int lvl_, int charge_) except + nogil 
        InternalCalibration_LockMass(InternalCalibration_LockMass &) except + nogil  # compiler
        double mz # m/z of the lock mass (incl. adducts)
        unsigned int ms_level # MS level where it occurs
        int charge # charge of the ion (to find isotopes)

