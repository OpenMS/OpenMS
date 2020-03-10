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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>

using namespace std;

namespace OpenMS
{

  SplinePackage::SplinePackage(std::vector<double> pos, std::vector<double> intensity) :
    spline_(pos, intensity)
  {
    if (!(pos.size() == intensity.size() && pos.size() > 1))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "m/z (or RT) and intensity vectors either not of the same size or too short.");
    }

    pos_min_ = pos.front();
    pos_max_ = pos.back();
    pos_step_width_ = (pos_max_ - pos_min_) / (pos.size() - 1);
  }

  SplinePackage::~SplinePackage()
  {
  }

  double SplinePackage::getPosMin() const
  {
    return pos_min_;
  }

  double SplinePackage::getPosMax() const
  {
    return pos_max_;
  }

  double SplinePackage::getPosStepWidth() const
  {
    return pos_step_width_;
  }

  bool SplinePackage::isInPackage(double pos) const
  {
    return pos >= pos_min_ && pos <= pos_max_;
  }

  double SplinePackage::eval(double pos) const
  {
    if (this->isInPackage(pos))
    {
      return max(0.0, spline_.eval(pos));
    }
    else
    {
      return 0;
    }
  }

}
