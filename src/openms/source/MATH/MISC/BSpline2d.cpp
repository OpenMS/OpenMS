// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Stephan Aiche, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MISC/BSpline2d.h>

#include <BSpline/BSplineBase.cpp>
#include <BSpline/BSpline.cpp>

namespace OpenMS
{

  BSpline2d::BSpline2d(const std::vector<double>& x, const std::vector<double>& y, double wave_length, BoundaryCondition boundary_condition, Size num_nodes)
  {
    OPENMS_PRECONDITION(x.size() == y.size(), "x and y vector passed to BSpline constructor must match.")
    spline_ = new BSpline<double>(&x[0], x.size(), &y[0], wave_length, boundary_condition, num_nodes);
  }

  BSpline2d::~BSpline2d()
  {
    delete spline_;
  }

  bool BSpline2d::solve (const std::vector<double>& y)
  {
    OPENMS_PRECONDITION(x.size() == y.size(), "y vector passed to BSpline solve must match size of x.")
    // pass vector as array
    return spline_->solve(&y[0]);
  }

  double BSpline2d::eval (double x)
  {
    return spline_->evaluate(x);
  }

  double BSpline2d::derivative (double x)
  {
    return spline_->slope(x);
  }

  double BSpline2d::coefficient (int n)
  {
    return spline_->coefficient(n);
  }

}
