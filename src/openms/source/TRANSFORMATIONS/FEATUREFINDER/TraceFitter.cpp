// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>

#include <unsupported/Eigen/NonLinearOptimization>

namespace OpenMS
{

  int TraceFitter::GenericFunctor::inputs() const
  {
    return m_inputs;
  }

  int TraceFitter::GenericFunctor::values() const
  {
    return m_values;
  }

  TraceFitter::GenericFunctor::GenericFunctor(int dimensions, int num_data_points) :
    m_inputs(dimensions), m_values(num_data_points)
  {
  }

  TraceFitter::GenericFunctor::~GenericFunctor()
  {
  }

  TraceFitter::TraceFitter() :
    DefaultParamHandler("TraceFitter")
  {
    defaults_.setValue("max_iteration", 500, "Maximum number of iterations used by the Levenberg-Marquardt algorithm.", ListUtils::create<String>("advanced"));
    defaults_.setValue("weighted", "false", "Weight mass traces according to their theoretical intensities.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("weighted", ListUtils::create<String>("true,false"));
    defaultsToParam_();
  }

  TraceFitter::TraceFitter(const TraceFitter& source) :
    DefaultParamHandler(source),
    max_iterations_(source.max_iterations_),
    weighted_(source.weighted_)
  {
    updateMembers_();
  }

  TraceFitter& TraceFitter::operator=(const TraceFitter& source)
  {
    DefaultParamHandler::operator=(source);
    max_iterations_ = source.max_iterations_;
    weighted_ = source.weighted_;
    updateMembers_();

    return *this;
  }

  TraceFitter::~TraceFitter()
  {
  }

  double TraceFitter::computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, Size k)
  {
    double rt = trace.peaks[k].first;

    return trace.theoretical_int * getValue(rt);
  }

  void TraceFitter::updateMembers_()
  {
    max_iterations_ = this->param_.getValue("max_iteration");
    weighted_ = this->param_.getValue("weighted") == "true";
  }

  void TraceFitter::optimize_(Eigen::VectorXd& x_init, GenericFunctor& functor)
  {
    //TODO: this function is copy&paste from LevMarqFitter1d.h. Make a generic wrapper for
    //LM optimization
    int data_count = functor.values();
    int num_params = functor.inputs();

    // LM always expects N>=p, cause Jacobian be rectangular M x N with M>=N
    if (data_count < num_params) throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-FinalSet", "Skipping feature, we always expects N>=p");


    Eigen::LevenbergMarquardt<GenericFunctor> lmSolver(functor);
    lmSolver.parameters.maxfev = max_iterations_;
    Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_init);

    //the states are poorly documented. after checking the source, we believe that
    //all states except NotStarted, Running and ImproperInputParameters are good
    //termination states.
    if (status <= Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-FinalSet", "Could not fit the gaussian to the data: Error " + String(status));
    }

    getOptimizedParameters_(x_init);
  }

} // namespace OpenMS
