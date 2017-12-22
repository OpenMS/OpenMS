// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_EGHTRACEFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_EGHTRACEFITTER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>

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

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDTRACEFITTERGAUSS_H
