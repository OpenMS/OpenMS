// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

namespace OpenMS
{

  GaussFilter::GaussFilter() :
    ProgressLogger(),
    DefaultParamHandler("GaussFilter"),
    coeffs_(),
    sigma_(0.1),
    spacing_(0.01)
  {
    //Parameter settings
    defaults_.setValue("gaussian_width", 0.2, "Use a gaussian filter width which has approximately the same width as your mass peaks (FWHM in m/z).");
    defaults_.setValue("ppm_tolerance", 10.0, "Gaussian width, depending on the m/z position.\nThe higher the value, the wider the peak and therefore the wider the gaussian.");
    defaults_.setValue("use_ppm_tolerance", "false", "If true, instead of the gaussian_width value, the ppm_tolerance is used. The gaussian is calculated in each step anew, so this is much slower.");
    defaults_.setValidStrings("use_ppm_tolerance", StringList::create("true,false"));
    defaultsToParam_();
  }

  GaussFilter::~GaussFilter()
  {
  }

  void GaussFilter::updateMembers_()
  {
    sigma_ = (DoubleReal)param_.getValue("gaussian_width") / 8.0;
    Size number_of_points_right = (Size)(ceil(4 * sigma_ / spacing_)) + 1;
    coeffs_.resize(number_of_points_right);
    coeffs_[0] = 1.0 / (sigma_ * sqrt(2.0 * Constants::PI));

    for (Size i = 1; i < number_of_points_right; i++)
    {
      coeffs_[i] = 1.0 / (sigma_ * sqrt(2.0 * Constants::PI)) * exp(-((i * spacing_) * (i * spacing_)) / (2 * sigma_ * sigma_));
    }
#ifdef DEBUG_FILTERING
    std::cout << "Coeffs: " << std::endl;
    for (Size i = 0; i < number_of_points_right; i++)
    {
      std::cout << i * spacing_ << ' ' << coeffs_[i] << std::endl;
    }
#endif
  }

}
