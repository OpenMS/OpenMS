// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  /**
    @brief Compute the area, background and shape metrics of a peak.

    The area computation is performed in integratePeak() and it supports
    integration by simple sum of the intensity, integration by Simpson's rule
    implementations for an odd number of unequally spaced points or integration
    by the trapezoid rule.

    The background computation is performed in estimateBackground() and it
    supports three different approaches to baseline correction, namely
    computing a rectangular shape under the peak based on the minimum value of
    the peak borders (vertical_division_min), a rectangular shape based on the
    maximum value of the beak borders (vertical_division_max) or a trapezoidal
    shape based on a straight line between the peak borders (base_to_base).

    Peak shape metrics are computed in calculatePeakShapeMetrics() and multiple
    metrics are supported.

    The containers supported by the methods are MSChromatogram and MSSpectrum.
  */
  class OPENMS_DLLAPI EmgGradientDescent :
    public DefaultParamHandler
  {
public:
    /// Constructor
    EmgGradientDescent();
    /// Destructor
    ~EmgGradientDescent() override = default;

    void getDefaultParameters(Param& params);

    /// To test private and protected methods
    friend class EmgGradientDescent_friend;

    /**
      @brief Fit the given peak (either MSChromatogram or MSSpectrum) to the EMG peak model

      The method is able to recapitulate the actual peak area of saturated or cutoff peaks.
      In addition, the method is able to fine tune the peak area of well acquired peaks.
      The output is a reconstruction of the input peak. Additional points are often added
      to produce a peak with similar intensities on boundaries' points.

      Metadata will be added to the output peak, containing the optimal parameters
      for the EMG peak model. This information will be found in a `FloatDataArray`
      of name "emg_parameters", with the parameters being saved in the following
      order (from index 0 to 3): amplitude `h`, mean `mu`, standard deviation `sigma`,
      exponent relaxation time `tau`.

      If `left_pos` and `right_pos` are passed, then only that part of the peak
      is taken into consideration.

      @note All optimal gradient descent parameters are currently hard coded to allow for a simplified user interface

      @note Cutoff peak: The intensities of the left and right baselines are not equal

      @note Saturated peak: The maximum intensity of the peak is lower than expected due to saturation of the detector

      Inspired by the results found in:
      Yuri Kalambet, Yuri Kozmin, Ksenia Mikhailova, Igor Nagaev, Pavel Tikhonov
      Reconstruction of chromatographic peaks using the exponentially modified Gaussian function

      @tparam PeakContainerT Either a MSChromatogram or a MSSpectrum
      @param[in] input_peak Input peak
      @param[out] output_peak Output peak
      @param[in] left_pos RT or MZ value of the first point of interest
      @param[in] right_pos RT or MZ value of the last point of interest
    */
    template <typename PeakContainerT>
    void fitEMGPeakModel(
      const PeakContainerT& input_peak,
      PeakContainerT& output_peak,
      const double left_pos = 0.0,
      const double right_pos = 0.0
    ) const;

    /**
      @brief The implementation of the gradient descent algorithm for the EMG peak model

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[out] best_h `h` (amplitude) parameter
      @param[out] best_mu `mu` (mean) parameter
      @param[out] best_sigma `sigma` (standard deviation) parameter
      @param[out] best_tau `tau` (exponent relaxation time) parameter

      @return The number of iterations necessary to reach the best values for the parameters
    */
    UInt estimateEmgParameters(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      double& best_h,
      double& best_mu,
      double& best_sigma,
      double& best_tau
    ) const;

    /**
      @brief Compute the EMG function on a set of points

      If class parameter `compute_additional_points` is `"true"`, the algorithm
      will detect which side of the peak is cutoff and add points to it.

      @param[in] xs Positions
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time
      @param[out] out_xs The output positions
      @param[out] out_ys The output intensities
    */
    void applyEstimatedParameters(
      const std::vector<double>& xs,
      const double h,
      const double mu,
      const double sigma,
      const double tau,
      std::vector<double>& out_xs,
      std::vector<double>& out_ys
    ) const;

protected:
    void updateMembers_() override;

    /**
      @brief Given a peak, extract a training set to be used with the gradient descent algorithm

      The algorithm tries to select only those points that can help in finding the optimal
      parameters with gradient descent. The decision of which points to skip is based on the
      derivatives between consecutive points.

      It first selects all those points whose intensity is below a certain value (`intensity_threshold`).
      Then, the derivatives of all the remaining points are computed. Based on the results,
      the algorithm selects those points that present a high enough derivative.
      Once a low value is found, the algorithm stops taking points from that side.
      It then repeats the same procedure on the other side of the peak.
      The goal is to limit the inclusion of saturated or spurious points near the
      peak apex during training.

      @throw Exception::SizeUnderflow if the input has less than 2 elements

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[out] TrX Extracted training set positions
      @param[out] TrY Extracted training set intensities
    */
    void extractTrainingSet(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      std::vector<double>& TrX,
      std::vector<double>& TrY
    ) const;

    /**
      @brief Compute the boundary for the mean (`mu`) parameter in gradient descent

      Together with the value returned by computeInitialMean(), this method
      decides the minimum and maximum value that `mu` can assume during iterations
      of the gradient descent algorithm.
      The value is based on the width of the peak.

      @param[in] xs Positions

      @return The maximum distance from the precomputed initial mean in the gradient descent algorithm
    */
    double computeMuMaxDistance(const std::vector<double>& xs) const;

    /**
      @brief Compute an estimation of the mean of a peak

      The method computes the middle point on different levels of intensity of the peak.
      The returned mean is the average of these middle points.

      @throw Exception::SizeUnderflow if the input is empty

      @param[in] xs Positions
      @param[in] ys Intensities

      @return The peak's estimated mean
    */
    double computeInitialMean(
      const std::vector<double>& xs,
      const std::vector<double>& ys
    ) const;

private:
    /**
      @brief Apply the iRprop+ algorithm for gradient descent

      Reference:
      Christian Igel and Michael Hüsken. Improving the Rprop Learning Algorithm.
      Second International Symposium on Neural Computation (NC 2000), pp. 115-121, ICSC Academic Press, 2000

      @param[in] prev_diff_E_param The cost of the partial derivative of E with
      respect to the given parameter, at the previous iteration of gradient descent
      @param[in,out] diff_E_param The cost of the partial derivative of E with
      respect to the given parameter, at the current iteration
      @param[in,out] param_lr The learning rate for the given parameter
      @param[in,out] param_update The amount to add/remove to/from `param`
      @param[in,out] param The parameter for which the algorithm tries speeding the convergence to a minimum
      @param[in] current_E The current cost E
      @param[in] previous_E The previous cost E
    */
    void iRpropPlus(
      const double prev_diff_E_param,
      double& diff_E_param,
      double& param_lr,
      double& param_update,
      double& param,
      const double current_E,
      const double previous_E
    ) const;

    /**
      @brief Compute the cost given by loss function E

      Needed by the gradient descent algorithm.
      The mean squared error is used as the loss function E.

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed cost
    */
    double Loss_function(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute the cost given by the partial derivative of the loss function E,
      with respect to `h` (the amplitude)

      Needed by the gradient descent algorithm.

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed cost
    */
    double E_wrt_h(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute the cost given by the partial derivative of the loss function E,
      with respect to `mu` (the mean)

      Needed by the gradient descent algorithm.

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed cost
    */
    double E_wrt_mu(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute the cost given by the partial derivative of the loss function E,
      with respect to `sigma` (the standard deviation)

      Needed by the gradient descent algorithm.

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed cost
    */
    double E_wrt_sigma(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute the cost given by the partial derivative of the loss function E,
      with respect to `tau` (the exponent relaxation time)

      Needed by the gradient descent algorithm.

      @param[in] xs Positions
      @param[in] ys Intensities
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed cost
    */
    double E_wrt_tau(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute EMG's z parameter

      The value of z decides which formula is to be used during EMG function computation.
      Z values in the following ranges will each use a different EMG formula to
      avoid numerical instability and potential numerical overflow:
      (-inf, 0), [0, 6.71e7], (6.71e7, +inf)

      Reference:
      Kalambet, Y.; Kozmin, Y.; Mikhailova, K.; Nagaev, I.; Tikhonov, P. (2011).
      "Reconstruction of chromatographic peaks using the exponentially modified
      Gaussian function". Journal of Chemometrics. 25 (7): 352.

      @param[in] x Position
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The computed parameter z
    */
    double compute_z(
      const double x,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /**
      @brief Compute the EMG function on a single point

      @param[in] x Position
      @param[in] h Amplitude
      @param[in] mu Mean
      @param[in] sigma Standard deviation
      @param[in] tau Exponent relaxation time

      @return The estimated intensity for the given input point
    */
    double emg_point(
      const double x,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const;

    /// Alias for OpenMS::Constants:PI
    const double PI = OpenMS::Constants::PI;

    /**
      Level of debug information to print to the terminal
      Valid values are: 0, 1, 2
      Higher values mean more information
    */
    UInt print_debug_;

    /// Maximum number of gradient descent iterations in `fitEMGPeakModel()`
    UInt max_gd_iter_;

    /**
      Whether additional points should be added when fitting EMG peak model,
      particularly useful with cutoff peaks
    */
    bool compute_additional_points_;
  };

  class EmgGradientDescent_friend
  {
public:
    EmgGradientDescent_friend() = default;
    ~EmgGradientDescent_friend() = default;

    double Loss_function(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const
    {
      return emg_gd_.Loss_function(xs, ys, h, mu, sigma, tau);
    }

    double computeMuMaxDistance(const std::vector<double>& xs) const
    {
      return emg_gd_.computeMuMaxDistance(xs);
    }

    void extractTrainingSet(
      const std::vector<double>& xs,
      const std::vector<double>& ys,
      std::vector<double>& TrX,
      std::vector<double>& TrY
    ) const
    {
      emg_gd_.extractTrainingSet(xs, ys, TrX, TrY);
    }

    double computeInitialMean(
      const std::vector<double>& xs,
      const std::vector<double>& ys
    ) const
    {
      return emg_gd_.computeInitialMean(xs, ys);
    }

    void iRpropPlus(
      const double prev_diff_E_param,
      double& diff_E_param,
      double& param_lr,
      double& param_update,
      double& param,
      const double current_E,
      const double previous_E
    ) const
    {
      emg_gd_.iRpropPlus(
        prev_diff_E_param, diff_E_param, param_lr,
        param_update, param, current_E, previous_E
      );
    }

    double compute_z(
      const double x,
      const double mu,
      const double sigma,
      const double tau
    ) const
    {
      return emg_gd_.compute_z(x, mu, sigma, tau);
    }

    void applyEstimatedParameters(
      const std::vector<double>& xs,
      const double h,
      const double mu,
      const double sigma,
      const double tau,
      std::vector<double>& out_xs,
      std::vector<double>& out_ys
    ) const
    {
      emg_gd_.applyEstimatedParameters(xs, h, mu, sigma, tau, out_xs, out_ys);
    }

    double emg_point(
      const double x,
      const double h,
      const double mu,
      const double sigma,
      const double tau
    ) const
    {
      return emg_gd_.emg_point(x, h, mu, sigma, tau);
    }

    EmgGradientDescent emg_gd_;
  };
}
