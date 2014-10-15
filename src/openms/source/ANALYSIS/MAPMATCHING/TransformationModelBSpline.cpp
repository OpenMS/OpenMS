// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>

using namespace std;

namespace OpenMS
{

  TransformationModelBSpline::TransformationModelBSpline(
    const TransformationModel::DataPoints& data, const Param& params)
  {
    // parameter handling/checking:
    params_ = params;
    Param defaults;
    getDefaultParameters(defaults);
    params_.setDefaults(defaults);

    if (data.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "'b_spline' model requires more data");
    }
    Size boundary_condition = params_.getValue("boundary_condition");

    BSpline2d::BoundaryCondition bound_cond = 
      static_cast<BSpline2d::BoundaryCondition>(boundary_condition);
    vector<double> x(data.size()), y(data.size());
    double xmin = data[0].first;
    double xmax = xmin;
    for (Size i = 0; i < data.size(); ++i)
    {
      x[i] = data[i].first;
      y[i] = data[i].second;
      if (x[i] < xmin) xmin = x[i];
      else if (x[i] > xmax) xmax = x[i];
    }
    double wavelength = params_.getValue("wavelength");
    if (wavelength > (xmax - xmin))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "B-spline 'wavelength' can't be larger than the data range (here: " + String(xmax - xmin) + ").", String(wavelength));
    }

    // since we can't initialize a BSpline2d object in the init list (no c'tor
    // that doesn't require preparation of data), we have to use a pointer:
    spline_ = new BSpline2d(x, y, wavelength, bound_cond, 
                            params_.getValue("num_nodes"));

    if (!spline_->ok())
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
                                   "TransformationModelBSpline", 
                                   "Unable to fit B-spline to data points.");
    }
  }


  TransformationModelBSpline::~TransformationModelBSpline()
  {
    delete spline_;
  }


  double TransformationModelBSpline::evaluate(double value)
  {
    return spline_->eval(value); // does this work for extrapolation?
  }


  void TransformationModelBSpline::getDefaultParameters(Param& params)
  {
    params.clear();
    params.setValue("wavelength", 0.0, "Determines the amount of smoothing by setting the number of nodes for the B-spline. The number is chosen so that the spline approximates a low-pass filter with this cutoff wavelength. The wavelength is given in the same units as the data; a higher value means more smoothing. '0' sets the number of nodes to twice the number of input points.");
    params.setMinFloat("wavelength", 0.0);
    params.setValue("num_nodes", 5, "Number of nodes for B-spline fitting. Overrides 'wavelength' if set (to two or greater). A lower value means more smoothing.");
    params.setMinInt("num_nodes", 0);
    params.setValue("boundary_condition", 2, "Boundary condition at endpoints: 0 (value zero), 1 (first derivative zero) or 2 (second derivative zero)");
    params.setMinInt("boundary_condition", 0);
    params.setMaxInt("boundary_condition", 2);
  }

} // namespace
