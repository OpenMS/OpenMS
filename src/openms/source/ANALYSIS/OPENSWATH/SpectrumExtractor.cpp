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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumExtractor.h>

namespace OpenMS
{
  SpectrumExtractor::SpectrumExtractor() :
    DefaultParamHandler("SpectrumExtractor")
  {
    getDefaultParameters(defaults_);

    // write defaults into Param object param_
    defaultsToParam_();
  }

  SpectrumExtractor::~SpectrumExtractor() {}

  void SpectrumExtractor::setRTWindow(const double& rt_window)
  {
    rt_window_ = rt_window;
  }

  double SpectrumExtractor::getRTWindow() const
  {
    return rt_window_;
  }

  void SpectrumExtractor::setMinScore(const double& min_score)
  {
    min_score_ = min_score;
  }

  double SpectrumExtractor::getMinScore() const
  {
    return min_score_;
  }

  void SpectrumExtractor::setMinForwardMatch(const double& min_forward_match)
  {
    min_forward_match_ = min_forward_match;
  }

  double SpectrumExtractor::getMinForwardMatch() const
  {
    return min_forward_match_;
  }

  void SpectrumExtractor::setMinReverseMatch(const double& min_reverse_match)
  {
    min_reverse_match_ = min_reverse_match;
  }

  double SpectrumExtractor::getMinReverseMatch() const
  {
    return min_reverse_match_;
  }

  void SpectrumExtractor::setMZTolerance(const double& mz_tolerance)
  {
    mz_tolerance_ = mz_tolerance;
  }

  double SpectrumExtractor::getMZTolerance() const
  {
    return mz_tolerance_;
  }

  void SpectrumExtractor::setMZToleranceUnits(const String& mz_tolerance_units)
  {
    mz_tolerance_units_ = mz_tolerance_units;
  }

  String SpectrumExtractor::getMZToleranceUnits() const
  {
    return mz_tolerance_units_;
  }

  void SpectrumExtractor::updateMembers_()
  {
    rt_window_ = (double)param_.getValue("rt_window");
    min_score_ = (double)param_.getValue("min_score");
    min_forward_match_ = (double)param_.getValue("min_forward_match");
    min_reverse_match_ = (double)param_.getValue("min_reverse_match");
    mz_tolerance_ = (double)param_.getValue("mz_tolerance");
    mz_tolerance_units_ = (String)param_.getValue("mz_tolerance_units");
  }

  void SpectrumExtractor::getDefaultParameters(Param& params)
  {
    params.clear();
    // TODO also set min and max values for these defaults (eg. params.setMinFloat(...))
    params.setValue("rt_window", 100, "Retention time window.");
    params.setValue("min_score", 0.7, "Minimum score.");
    params.setValue("min_forward_match", 0.7, "Minimum forward match.");
    params.setValue("min_reverse_match", 0.7, "Minimum reverse match.");
    params.setValue("mz_tolerance", 0.7, "Mass to Charge tolerance.");
    params.setValue("mz_tolerance_units", "ppm", "Mass to Charge tolerance units.");
    params.setValidStrings("mz_tolerance_units", ListUtils::create<String>("ppm,Da"));
  }
}
