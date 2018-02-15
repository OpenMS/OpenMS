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
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>

namespace OpenMS
{
  PeakShape::PeakShape() 
    : height(0),
      mz_position(0),
      left_width(0),
      right_width(0),
      area(0),
      r_value(0),
      signal_to_noise(0.),
      type(UNDEFINED),
      left_iterator_set_(false),
      right_iterator_set_(false)
  {
    left_endpoint_ = exp_.end();
    right_endpoint_ = exp_.end();
  }

  PeakShape::PeakShape(double height_arg, double mz_position_arg, double left_width_arg, 
      double right_width_arg, double area_arg, PeakIterator left_arg, PeakIterator right_arg, Type type_arg) :
    height(height_arg),
    mz_position(mz_position_arg),
    left_width(left_width_arg),
    right_width(right_width_arg),
    area(area_arg),
    r_value(0),
    signal_to_noise(0),
    type(type_arg),
    left_endpoint_(left_arg),
    right_endpoint_(right_arg),
    left_iterator_set_(true),
    right_iterator_set_(true)
  {
  }

  PeakShape::PeakShape(double height_arg, double mz_position_arg, double left_width_arg, double right_width_arg, double area_arg, Type type_arg) :
    height(height_arg),
    mz_position(mz_position_arg),
    left_width(left_width_arg),
    right_width(right_width_arg),
    area(area_arg),
    r_value(0),
    signal_to_noise(0),
    type(type_arg),
    left_iterator_set_(false),
    right_iterator_set_(false)
  {
    left_endpoint_ = exp_.end();
    right_endpoint_ = exp_.end();
  }

  PeakShape::PeakShape(const PeakShape & rhs) :
    height(rhs.height),
    mz_position(rhs.mz_position),
    left_width(rhs.left_width),
    right_width(rhs.right_width),
    area(rhs.area),
    r_value(rhs.r_value),
    signal_to_noise(rhs.signal_to_noise),
    type(rhs.type),
    left_iterator_set_(rhs.left_iterator_set_),
    right_iterator_set_(rhs.right_iterator_set_)
  {
    if (rhs.iteratorsSet())
    {
      left_endpoint_ = rhs.getLeftEndpoint();
      right_endpoint_ = rhs.getRightEndpoint();
      left_iterator_set_ = true;
      right_iterator_set_ = true;
    }
    else
    {
      left_endpoint_ = exp_.end();
      right_endpoint_ = exp_.end();
    }

  }

  PeakShape::~PeakShape()
  {
  }

  PeakShape & PeakShape::operator=(const PeakShape & rhs)
  {
    // handle self assignment
    if (this == &rhs)
      return *this;

    height = rhs.height;
    mz_position = rhs.mz_position;
    left_width = rhs.left_width;
    right_width = rhs.right_width;
    area = rhs.area;
    type = rhs.type;
    signal_to_noise = rhs.signal_to_noise;
    if (rhs.iteratorsSet())
    {
      left_endpoint_ = rhs.getLeftEndpoint();
      right_endpoint_ = rhs.getRightEndpoint();
      left_iterator_set_ = true;
      right_iterator_set_ = true;
    }
    else
    {
      left_endpoint_ = exp_.end();
      right_endpoint_ = exp_.end();
    }
    r_value = rhs.r_value;

    return *this;
  }

  bool PeakShape::operator==(const PeakShape & rhs) const
  {
    return height == rhs.height &&
           mz_position == rhs.mz_position &&
           left_width == rhs.left_width &&
           right_width == rhs.right_width &&
           area == rhs.area &&
           type == rhs.type &&
           signal_to_noise == rhs.signal_to_noise &&
           r_value == rhs.r_value;
  }

  bool PeakShape::operator!=(const PeakShape & rhs) const
  {
    return !(*this==rhs);
  }

  double PeakShape::operator()(double x) const
  {
    double value;

    switch (type)
    {
    case LORENTZ_PEAK:
      if (x <= mz_position)
        // see equation 8.7 on p. 72 in Dissertation of Eva Lange
        value = height / (1. + pow(left_width * (x - mz_position), 2));
      else
        value = height / (1. + pow(right_width * (x - mz_position), 2));
      break;

    case SECH_PEAK:
      if (x <= mz_position)
        // see equation 8.8 on p. 72 in Dissertation of Eva Lange
        value = height / pow(cosh(left_width * (x - mz_position)), 2);
      else
        value = height / pow(cosh(right_width * (x - mz_position)), 2);
      break;

    default:
      value = -1.;
      break;
    }

    return value;
  }

  double PeakShape::getFWHM() const
  {
    // static double m = log1p(sqrt(2.0));
    const double m = 0.88137358701954305; // log(1+sqrt(2.0))

    double fwhm = 0;
    if (right_width == 0. || left_width == 0.)
    {
      return -1.;
    }

    switch (type)
    {
      case LORENTZ_PEAK:
        // see equation 8.20 on p. 75 in Dissertation of Eva Lange
        fwhm = 1 / right_width;
        fwhm += 1 / left_width;
        break;

      case SECH_PEAK:
        // see equation 8.22 on p. 75 in Dissertation of Eva Lange
        fwhm = m / left_width;
        fwhm += m / right_width;
        break;

      default:
        fwhm = -1.;
        break;
    }
    return fwhm;
  }

  double PeakShape::getSymmetricMeasure() const
  {
    double value;

    if (left_width < right_width)
      value = left_width / right_width;
    else
      value = right_width / left_width;

    return value;
  }

  bool PeakShape::iteratorsSet() const
  {
    if (left_iterator_set_ && right_iterator_set_)
      return true;
    else
      return false;
  }

  PeakShape::PeakIterator PeakShape::getLeftEndpoint() const
  {
    return left_endpoint_;
  }

  void PeakShape::setLeftEndpoint(PeakShape::PeakIterator left_endpoint)
  {
    left_endpoint_ = left_endpoint;
    left_iterator_set_ = true;
  }

  PeakShape::PeakIterator PeakShape::getRightEndpoint() const
  {
    return right_endpoint_;
  }

  void PeakShape::setRightEndpoint(PeakShape::PeakIterator right_endpoint)
  {
    right_endpoint_ = right_endpoint;
    right_iterator_set_ = true;
  }

}
