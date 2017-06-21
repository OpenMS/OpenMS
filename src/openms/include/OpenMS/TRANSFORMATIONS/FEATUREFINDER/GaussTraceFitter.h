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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Marc Sturm$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_GAUSSTRACEFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_GAUSSTRACEFITTER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>

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

    virtual ~GaussTraceFitter();

    // override important methods
    void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces);

    double getLowerRTBound() const;

    double getUpperRTBound() const;

    double getHeight() const;

    double getCenter() const;

    double getFWHM() const;

    /**
     * @brief Returns the sigma of the fitted gaussian model
     */
    double getSigma() const;

    bool checkMaximalRTSpan(const double max_rt_span);

    bool checkMinimalRTSpan(const std::pair<double, double>& rt_bounds, const double min_rt_span);

    double getValue(double rt) const;

    double getArea();

    String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, const char function_name, const double baseline, const double rt_shift);

protected:
    double sigma_;
    double x0_;
    double height_;
    double region_rt_span_;

    static const Size NUM_PARAMS_;

    void getOptimizedParameters_(const Eigen::VectorXd& x_init);

    class GaussTraceFunctor :
      public TraceFitter::GenericFunctor
    {
public:
      GaussTraceFunctor(int dimensions,
                        const TraceFitter::ModelData* data);

      int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec);

      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J);
protected:
      const TraceFitter::ModelData* m_data;
    };

    void setInitialParameters_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces);

    virtual void updateMembers_();

  };

} // namespace OpenMS

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDTRACEFITTERGAUSS_H
