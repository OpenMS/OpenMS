// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/TraceFitter.h>

namespace OpenMS
{

  /**
   * @brief A RT Profile fitter using an Exponential Gaussian Hybrid background model
   *
   * Lan K, Jorgenson JW.
   * <b>A hybrid of exponential and gaussian functions as a simple model of asymmetric chromatographic peaks.</b>
   * <em>Journal of Chromatography A.</em> 915 (1-2)p. 1-13.
   * Available at: http://linkinghub.elsevier.com/retrieve/pii/S0021967301005945
   *
   * @htmlinclude OpenMS_EGHTraceFitter.parameters
   *
   * @experimental Needs further testing on real data. Note that the tests are currently also focused on testing the EGH as replacement for the gaussian.
   */
  class OPENMS_DLLAPI EGHTraceFitter :
    public TraceFitter
  {
public:
    /** Functor for LM Optimization */
    class EGHTraceFunctor :
      public TraceFitter::GenericFunctor
    {
public:

      EGHTraceFunctor(int dimensions,
                      const TraceFitter::ModelData* data);

      ~EGHTraceFunctor() override;

      int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) override;

      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) override;

protected:
      const TraceFitter::ModelData* m_data;
    };

    EGHTraceFitter();

    EGHTraceFitter(const EGHTraceFitter& other);

    EGHTraceFitter& operator=(const EGHTraceFitter& source);

    ~EGHTraceFitter() override;

    // override important methods
    void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces) override;

    double getLowerRTBound() const override;

    double getTau() const;

    double getUpperRTBound() const override;

    double getHeight() const override;

    double getSigma() const;

    double getCenter() const override;

    bool checkMaximalRTSpan(const double max_rt_span) override;

    bool checkMinimalRTSpan(const std::pair<double, double>& rt_bounds, const double min_rt_span) override;

    double getValue(double rt) const override;

    double getArea() override;

    double getFWHM() const override;

    String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, const char function_name, const double baseline, const double rt_shift) override;

protected:
    double apex_rt_;
    double height_;

    double sigma_;
    double tau_;

    std::pair<double, double> sigma_5_bound_;

    double region_rt_span_;

    /// Coefficients to calculate the proportionality factor for the peak area
    static const double EPSILON_COEFS_[];

    static const Size NUM_PARAMS_;

    /**
     * @brief Return an ordered pair of the positions where the EGH reaches a height of alpha * height of the EGH
     *
     * @param alpha The alpha at which the boundaries should be computed
     */
    std::pair<double, double> getAlphaBoundaries_(const double alpha) const;

    void getOptimizedParameters_(const Eigen::VectorXd& x_init) override;

    void setInitialParameters_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces);

    void updateMembers_() override;
  };

} // namespace OpenMS

