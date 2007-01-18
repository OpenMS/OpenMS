// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>

namespace OpenMS
{
  PeakShape::PeakShape(double height_,
                       double mz_position_,
                       double left_width_,
                       double right_width_,
                       double area_,
                       PeakShapeType::Enum type_)
      : height(height_),
      mz_position(mz_position_),
      left_width(left_width_),
      right_width(right_width_),
      area(area_),
      r_value(0),
			signal_to_noise(0),
			type(type_)
  {}

  PeakShape::PeakShape(const PeakShape& peakshape)
      : height(peakshape.height),
      mz_position(peakshape.mz_position),
      left_width(peakshape.left_width),
      right_width(peakshape.right_width),
      area(peakshape.area),
      r_value(peakshape.r_value),
      signal_to_noise(peakshape.signal_to_noise),
      type(peakshape.type)
  {}

  PeakShape& PeakShape::operator = (const PeakShape& pf)
  {
    // handle self assignment
    if (this == &pf) return *this;

    height=pf.height;
    mz_position=pf.mz_position;
    left_width=pf.left_width;
    right_width=pf.right_width;
    area=pf.area;
    type=pf.type;
    signal_to_noise = pf.signal_to_noise;
    r_value=pf.r_value;

    return *this;
  }

  double PeakShape::operator () (const double x) const
  {
    double value;

    switch (type)
    {
    case PeakShapeType::LORENTZ_PEAK:
      if (x<=mz_position)
        value = height/(1.+pow(left_width*(x - mz_position), 2));
      else
        value = height/(1.+pow(right_width*(x - mz_position), 2));
      break;

    case PeakShapeType::SECH_PEAK:
      if (x<=mz_position)
        value = height/pow(cosh(left_width*(x-mz_position)), 2);
      else
        value = height/pow(cosh(right_width*(x-mz_position)), 2);
      break;

    default:
      value = -1.;
      break;
    }

    return value;
  }

  double PeakShape::getFWHM() const
  {
    double fwhm=0;

    switch (type)
    {
    case PeakShapeType::LORENTZ_PEAK:
      {
        fwhm = 1/right_width;
        fwhm += 1/left_width;
      }
      break;

    case PeakShapeType::SECH_PEAK:
      {
        double m = log(sqrt(2.0)+1);
        fwhm = m/left_width;
        fwhm += m/right_width;
      }
      break;

    default:
      {
        fwhm = -1.;
      }
      break;
    }
    return fwhm;
  }

  double PeakShape::getSymmetricMeasure() const
  {
    double value;

    if (left_width < right_width)
      value = left_width/right_width;
    else
      value =	right_width/left_width;

    return value;
  }
}
