// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/CalibrationData.h>
#include <OpenMS/ML/RANSAC/RANSAC.h>

#include <vector>

namespace OpenMS
{
    
  /**

    @brief Create and apply models of a mass recalibration function.

    The input is a list of calibration points (ideally spanning a wide m/z range to prevent extrapolation when applying to model).
    
    Models (LINEAR, LINEAR_WEIGHTED, QUADRATIC, QUADRATIC_WEIGHTED) can be trained using CalData points (or a subset of them).
    Calibration points can have different retention time points, and a model should be build such that it captures
    the local (in time) decalibration of the instrument, i.e. choose appropriate time windows along RT to calibrate the
    spectra in this RT region.
    From the available calibrant data, a model is build. Later, any uncalibrated m/z value can be fed to the model, to obtain
    a calibrated m/z.

    The input domain can either be absolute mass differences in [Th], or relative differences in [ppm].
    The models are build based on this input.

    Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models.

  */
  class OPENMS_DLLAPI MZTrafoModel
  {
  
  private:
    std::vector<double> coeff_;  ///< Model coefficients (for both linear and quadratic models), estimated from the data
    bool use_ppm_; ///< during training, model is build on absolute or relative(ppm) predictions. predict(), i.e. applying the model, requires this information too
    double rt_; ///< retention time associated to the model (i.e. where the calibrant data was taken from)

    static Math::RANSACParam* ransac_params_; ///< global pointer, init to NULL at startup; set class-global RANSAC params
    static int ransac_seed_; ///< seed used for all RANSAC invocations
    static double limit_offset_; ///< acceptable boundary for the estimated offset; if estimated offset is larger (absolute) the model does not validate (isValidModel())
    static double limit_scale_; ///< acceptable boundary for the estimated scale; if estimated scale is larger (absolute) the model does not validate (isValidModel())
    static double limit_power_; ///< acceptable boundary for the estimated power; if estimated power is larger (absolute) the model does not validate (isValidModel())

  public:

    /**
      @brief Default constructor
    */
    MZTrafoModel();

    /**
      @brief Default constructor

      If you have external coefficients, use this constructor and the setCoefficients() method to
      build a 'manual' model.
      Afterwards, use applyTransformation() or predict() to calibrate your data.
      If you call train(), the ppm-setting will be overwritten, depending on the type of training data.

      @param ppm_model Are the coefficients derived from ppm calibration data, or from absolute deltas?
    */
    MZTrafoModel(bool ppm_model);

    enum MODELTYPE { LINEAR, LINEAR_WEIGHTED, QUADRATIC, QUADRATIC_WEIGHTED, SIZE_OF_MODELTYPE };
    static const std::string names_of_modeltype[]; ///< strings corresponding to enum MODELTYPE
    /**
      @brief Convert string to enum

      Returns 'SIZE_OF_MODELTYPE' if string is unknown.
      @param name A string from names_of_modeltype[].
      @return The corresponding enum value.
    */
    static MODELTYPE nameToEnum(const std::string& name);
    /**
      @brief Convert enum to string
        
      @param mt The enum value
      @return Stringified version
    */
    static const std::string& enumToName(MODELTYPE mt); 


    /**
      @brief Set the global (program wide) parameters for RANSAC.

      This is not done via member, to keep a small memory footprint since hundreds of
      MZTrafoModels are expected to be build at the same time and the RANSAC params
      should be identical for all of them.
      
      @param p RANSAC params
    */
    static void setRANSACParams(const Math::RANSACParam& p);

    /**
      @brief Set RANSAC seed
    */
    static void setRANSACSeed(int seed);

    /**
      @brief Set coefficient boundaries for which the model coefficient must not exceed to be considered a valid model

      Use std::numeric_limits<double>::max() for no limit (default).
      If isValidModel() is called these limits are checked.
      Negative input run through fabs() to get positive values (since comparison is done in absolute terms).

    */
    static void setCoefficientLimits(double offset, double scale, double power);
    /**
      @brief Predicate to decide if the model has valid parameters, i.e. coefficients.

      If the model coefficients are empty, no model was trained yet (or unsuccessful),
      causing a return value of 'false'.

      Also, if the model has coefficients, we check if they are within the
      acceptable boundaries (if boundaries were given via setCoeffientLimits()).

    */
    static bool isValidModel(const MZTrafoModel& trafo);

    /**
      @brief Does the model have coefficients (i.e. was trained successfully).

      Having coefficients does not mean its valid (see isValidModel(); since coeffs might be too large).

    */
    bool isTrained() const;

    /**
      @brief Get RT associated with the model (training region)
    */
    double getRT() const;

    /**
      @brief Apply the model to an uncalibrated m/z value.

      Make sure the model was trained (train()) and is valid (isValidModel()) before calling this function!

      Applies the function y = intercept + slope*mz + power*mz^2
      and returns y.

      @param mz The uncalibrated m/z value
      @return The calibrated m/z value

    */
    double predict(double mz) const;
    
    /**
      @brief Binary search for the model nearest to a specific RT

      @param tms Vector of models, sorted by RT
      @param rt The target retention time
      @return Returns the index into 'tms' with the closest RT.

      @note Make sure the vector is sorted with respect to RT! Otherwise the result is undefined.

      @exception Exception::Precondition is thrown if the vector is empty (not only in debug mode)
    */
    static Size findNearest(const std::vector<MZTrafoModel>& tms, double rt);


    /// Comparator by position. As this class has dimension 1, this is basically an alias for MZLess.
    struct RTLess
    {
      inline bool operator()(const double& left, const MZTrafoModel& right) const
      {
        return left < right.rt_;
      }
      inline bool operator()(const MZTrafoModel& left, const double& right) const
      {
        return left.rt_ < right;
      }
      inline bool operator()(const MZTrafoModel& left, const MZTrafoModel& right) const
      {
        return left.rt_ < right.rt_;
      }
    };
    /**
      @brief Train a model using calibrant data

      If the CalibrationData was created using peak groups (usually corresponding to mass traces),
      the median for each group is used as a group representative. This
      is more robust, and reduces the number of data points drastically, i.e. one value per group.

      Internally, these steps take place:
       - apply RT filter
       - [compute median per group] (only if groups were given in 'cd')
       - set Model's rt position
       - call train() (see overloaded method)

      @param cd List of calibrants
      @param md Type of model (linear, quadratic, ...)
      @param use_RANSAC Remove outliers before computing the model?
      @param rt_left Filter 'cd' by RT; all calibrants with RT < 'rt_left' are removed
      @param rt_right Filter 'cd' by RT; all calibrants with RT > 'rt_right' are removed
      @return True if model was build, false otherwise

    */
    bool train(const CalibrationData& cd, MODELTYPE md, bool use_RANSAC, 
               double rt_left = -std::numeric_limits<double>::max(), 
               double rt_right = std::numeric_limits<double>::max()
               );

    /**
      @brief Train a model using calibrant data

      Given theoretical and observed mass values (and corresponding weights),
      a model (linear, quadratic, ...) is build.
      Outlier removal is applied before.
      The 'obs_mz' can be either given as absolute masses in [Th] or relative deviations in [ppm].
      The MZTrafoModel must be constructed accordingly (see constructor). This has no influence on the model building itself, but
      rather on how 'predict()' works internally.

      Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models.

      Internally, these steps take place:
       - [apply RANSAC] (depending on 'use_RANSAC')
       - build model and store its parameters internally

      @param error_mz Observed Mass error (in ppm or Th)
      @param theo_mz Theoretical m/z values, corresponding to 'error_mz'
      @param weights For weighted models only: weight of calibrants; ignored otherwise
      @param md Type of model (linear, quadratic, ...)
      @param use_RANSAC Remove outliers before computing the model?
      @return True if model was build, false otherwise

    */
    bool train(std::vector<double> error_mz,
               std::vector<double> theo_mz,
               std::vector<double> weights,
               MODELTYPE md,
               bool use_RANSAC);

    /**
      @brief Get model coefficients.

      Parameters will be filled with internal model parameters.
      The model must be trained before; Exception is thrown otherwise!

      @param intercept The intercept
      @param slope The slope
      @param power The coefficient for x*x (will be 0 for linear models)
      @throw Exception::Precondition if model is not trained yet
    */
    void getCoefficients(double& intercept, double& slope, double& power);

    /**
      @brief Copy model coefficients from another model.
    */
    void setCoefficients(const MZTrafoModel& rhs);

    /**
      @brief Manually set model coefficients

      Can be used instead of train(), so manually set coefficients.
      It must be exactly three values. If you want a linear model, set 'power' to zero.
      If you want a constant model, set slope to zero in addition.

      @param intercept The offset
      @param slope The slope
      @param power The x*x coefficient (for quadratic models)
    */
    void setCoefficients(double intercept, double slope, double power);

    /**
      @brief String representation of the model parameters.

      Empty if model is not trained.

    */
    String toString() const;

  }; // MZTrafoModel

} // namespace OpenMS

