// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Marc Sturm$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/TraceFitter.h>

namespace OpenMS
{

  /**
   * @brief Fitter for RT profiles using a Gaussian background model
   *
   * @htmlinclude OpenMS_GaussTraceFitter.parameters
   *
   * @todo More docu
   */
  class OPENMS_DLLAPI GaussTraceFitter :
    public TraceFitter
  {
public:
    GaussTraceFitter();

    GaussTraceFitter(const GaussTraceFitter& other);

    GaussTraceFitter& operator=(const GaussTraceFitter& source);

    ~GaussTraceFitter() override;

    // override important methods
    void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces) override;

    double getLowerRTBound() const override;

    double getUpperRTBound() const override;

    double getHeight() const override;

    double getCenter() const override;

    double getFWHM() const override;

    /**
     * @brief Returns the sigma of the fitted gaussian model
     */
    double getSigma() const;

    bool checkMaximalRTSpan(const double max_rt_span) override;

    bool checkMinimalRTSpan(const std::pair<double, double>& rt_bounds, const double min_rt_span) override;

    double getValue(double rt) const override;

    double getArea() override;

    String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, const char function_name, const double baseline, const double rt_shift) override;

protected:
    double sigma_;
    double x0_;
    double height_;
    double region_rt_span_;

    static const Size NUM_PARAMS_;

    void getOptimizedParameters_(const Eigen::VectorXd& s) override;

    class GaussTraceFunctor :
      public TraceFitter::GenericFunctor
    {
public:
      GaussTraceFunctor(int dimensions,
                        const TraceFitter::ModelData* data);

      int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) override;

      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) override;
protected:
      const TraceFitter::ModelData* m_data;
    };

    void setInitialParameters_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces);

    void updateMembers_() override;

  };

} // namespace OpenMS

