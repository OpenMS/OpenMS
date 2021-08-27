from Types cimport *
from MZTrafoModel cimport *
from RANSAC cimport *
from String cimport *
from CalibrationData cimport *

cdef extern from "<OpenMS/FILTERING/CALIBRATION/MZTrafoModel.h>" namespace "OpenMS":

    cdef cppclass MZTrafoModel:
        # wrap-doc:
                #   Create and apply models of a mass recalibration function
                #   -----
                #   The input is a list of calibration points (ideally spanning a wide m/z range to prevent extrapolation when applying to model)
                #   -----
                #   Models (LINEAR, LINEAR_WEIGHTED, QUADRATIC, QUADRATIC_WEIGHTED) can be trained using CalData points (or a subset of them)
                #   Calibration points can have different retention time points, and a model should be build such that it captures
                #   the local (in time) decalibration of the instrument, i.e. choose appropriate time windows along RT to calibrate the
                #   spectra in this RT region
                #   From the available calibrant data, a model is build. Later, any uncalibrated m/z value can be fed to the model, to obtain
                #   a calibrated m/z
                #   -----
                #   The input domain can either be absolute mass differences in [Th], or relative differences in [ppm]
                #   The models are build based on this input
                #   -----
                #   Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models

        MZTrafoModel()  nogil except +
        MZTrafoModel(MZTrafoModel &) nogil except + # compiler
        MZTrafoModel(bool) nogil except +
        bool isTrained() nogil except + # wrap-doc:Returns true if the model have coefficients (i.e. was trained successfully)
        double getRT() nogil except + # wrap-doc:Get RT associated with the model (training region)
        double predict(double mz) nogil except +
            # wrap-doc:
                #   Apply the model to an uncalibrated m/z value
                #   -----
                #   Make sure the model was trained (train()) and is valid (isValidModel()) before calling this function!
                #   -----
                #   Applies the function y = intercept + slope*mz + power*mz^2
                #   and returns y
                #   -----
                #   :param mz: The uncalibrated m/z value
                #   :returns The calibrated m/z value

        bool train(CalibrationData cd, MZTrafoModel_MODELTYPE md, bool use_RANSAC, double rt_left, double rt_right) nogil except +
            # wrap-doc:
                #   Train a model using calibrant data
                #   -----
                #   If the CalibrationData was created using peak groups (usually corresponding to mass traces),
                #   the median for each group is used as a group representative. This
                #   is more robust, and reduces the number of data points drastically, i.e. one value per group
                #   -----
                #   Internally, these steps take place:
                #   - apply RT filter
                #   - [compute median per group] (only if groups were given in 'cd')
                #   - set Model's rt position
                #   - call train() (see overloaded method)
                #   -----
                #   :param cd: List of calibrants
                #   :param md: Type of model (linear, quadratic, ...)
                #   :param use_RANSAC: Remove outliers before computing the model?
                #   :param rt_left: Filter 'cd' by RT; all calibrants with RT < 'rt_left' are removed
                #   :param rt_right: Filter 'cd' by RT; all calibrants with RT > 'rt_right' are removed
                #   :returns: True if model was build, false otherwise
                
        bool train(libcpp_vector[double] error_mz, libcpp_vector[double] theo_mz, libcpp_vector[double] weights, MZTrafoModel_MODELTYPE md, bool use_RANSAC) nogil except +
            # wrap-doc:
                #   Train a model using calibrant data
                #   -----
                #   Given theoretical and observed mass values (and corresponding weights),
                #   a model (linear, quadratic, ...) is build
                #   Outlier removal is applied before
                #   The 'obs_mz' can be either given as absolute masses in [Th] or relative deviations in [ppm]
                #   The MZTrafoModel must be constructed accordingly (see constructor). This has no influence on the model building itself, but
                #   rather on how 'predict()' works internally
                #   -----
                #   Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models
                #   -----
                #   Internally, these steps take place:
                #   - [apply RANSAC] (depending on 'use_RANSAC')
                #   - build model and store its parameters internally
                #   -----
                #   :param error_mz: Observed Mass error (in ppm or Th)
                #   :param theo_mz: Theoretical m/z values, corresponding to 'error_mz'
                #   :param weights: For weighted models only: weight of calibrants; ignored otherwise
                #   :param md: Type of model (linear, quadratic, ...)
                #   :param use_RANSAC: Remove outliers before computing the model?
                #   :returns: True if model was build, false otherwise

        void getCoefficients(double& intercept, double& slope, double& power) nogil except +
            # wrap-doc:
                #   Get model coefficients
                #   -----
                #   Parameters will be filled with internal model parameters
                #   The model must be trained before; Exception is thrown otherwise!
                #   -----
                #   :param intercept: The intercept
                #   :param slope: The slope
                #   :param power: The coefficient for x*x (will be 0 for linear models)

        void setCoefficients(MZTrafoModel) nogil except + # wrap-doc:Copy model coefficients from another model
        void setCoefficients(double, double, double) nogil except +
            # wrap-doc:
                #   Manually set model coefficients
                #   -----
                #   Can be used instead of train(), so manually set coefficients
                #   It must be exactly three values. If you want a linear model, set 'power' to zero
                #   If you want a constant model, set slope to zero in addition
                #   -----
                #   :param intercept: The offset
                #   :param slope: The slope
                #   :param power: The x*x coefficient (for quadratic models)

        String toString() nogil except +
        
cdef extern from "<OpenMS/FILTERING/CALIBRATION/MZTrafoModel.h>" namespace "OpenMS::MZTrafoModel":

    cdef enum MZTrafoModel_MODELTYPE "OpenMS::MZTrafoModel::MODELTYPE":
        LINEAR
        LINEAR_WEIGHTED
        QUADRATIC
        QUADRATIC_WEIGHTED
        SIZE_OF_MODELTYPE
    
    # static members
    # libcpp_string names_of_modeltype[] nogil except +
    MZTrafoModel_MODELTYPE nameToEnum(libcpp_string name) nogil except + # wrap-attach:MZTrafoModel
    libcpp_string enumToName(MZTrafoModel_MODELTYPE mt) nogil except + # wrap-attach:MZTrafoModel
    void setRANSACParams(RANSACParam p) nogil except + # wrap-attach:MZTrafoModel
    void setCoefficientLimits(double offset, double scale, double power) nogil except + # wrap-attach:MZTrafoModel
    bool isValidModel(MZTrafoModel& trafo) nogil except + # wrap-attach:MZTrafoModel
    Size findNearest(libcpp_vector[MZTrafoModel]& tms, double rt) nogil except + # wrap-attach:MZTrafoModel
    
