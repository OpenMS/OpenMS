// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>

namespace OpenMS
{
  PeakShape::PeakShape(DoubleReal height, DoubleReal mz_position, DoubleReal left_width, DoubleReal right_width, DoubleReal area, PeakIterator left, PeakIterator right, Type type)
		: height(height),
      mz_position(mz_position),
      left_width(left_width),
      right_width(right_width),
      area(area),
      r_value(0),
			signal_to_noise(0),
			type(type),
			left_endpoint_(left),
			right_endpoint_(right),
			left_iterator_set_(true),
			right_iterator_set_(true)
  {
  }

PeakShape::PeakShape(DoubleReal height, DoubleReal mz_position, DoubleReal left_width, DoubleReal right_width, DoubleReal area, Type type)
		: height(height),
      mz_position(mz_position),
      left_width(left_width),
      right_width(right_width),
      area(area),
      r_value(0),
			signal_to_noise(0),
			type(type),
			left_iterator_set_(false),
			right_iterator_set_(false)
  {
		left_endpoint_ = exp_.end();
		right_endpoint_ = exp_.end();
  }

	
  PeakShape::PeakShape(const PeakShape& rhs)
		: height(rhs.height),
      mz_position(rhs.mz_position),
      left_width(rhs.left_width),
      right_width(rhs.right_width),
      area(rhs.area),
      r_value(rhs.r_value),
      signal_to_noise(rhs.signal_to_noise),
      type(rhs.type)
  {
		if(rhs.iteratorsSet())
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

  PeakShape& PeakShape::operator = (const PeakShape& rhs)
  {
    // handle self assignment
    if (this == &rhs) return *this;

    height=rhs.height;
    mz_position=rhs.mz_position;
    left_width=rhs.left_width;
    right_width=rhs.right_width;
    area=rhs.area;
    type=rhs.type;
    signal_to_noise = rhs.signal_to_noise;
		if(rhs.iteratorsSet())
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
    r_value=rhs.r_value;

    return *this;
  }

	bool PeakShape::operator==(const PeakShape& rhs) const
	{
		return
			height==rhs.height && 
			mz_position==rhs.mz_position && 
			left_width==rhs.left_width && 
			right_width==rhs.right_width && 
			area==rhs.area && 
			type==rhs.type && 
			signal_to_noise == rhs.signal_to_noise && 
			r_value==rhs.r_value;
	}

	bool PeakShape::operator!=(const PeakShape& rhs) const
	{
		return
			height!=rhs.height ||
			mz_position!=rhs.mz_position || 
			left_width!=rhs.left_width || 
			right_width!=rhs.right_width || 
			area!=rhs.area || 
			type!=rhs.type || 
			signal_to_noise != rhs.signal_to_noise || 
			r_value!=rhs.r_value;
	}

  DoubleReal PeakShape::operator () (DoubleReal x) const
  {
    DoubleReal value;

    switch (type)
    {
    case LORENTZ_PEAK:
      if (x<=mz_position)
        value = height/(1.+pow(left_width*(x - mz_position), 2));
      else
        value = height/(1.+pow(right_width*(x - mz_position), 2));
      break;

    case SECH_PEAK:
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

  DoubleReal PeakShape::getFWHM() const
  {
    DoubleReal fwhm=0;
		if(right_width == 0. || left_width == 0.)
		{
				return -1.;
		}
    
    switch (type)
    {
    case LORENTZ_PEAK:
      {
		    fwhm = 1/right_width;
        fwhm += 1/left_width;
      }
      break;

    case SECH_PEAK:
      {
        DoubleReal m = log(sqrt(2.0)+1);
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

  DoubleReal PeakShape::getSymmetricMeasure() const
  {
    DoubleReal value;

    if (left_width < right_width)
      value = left_width/right_width;
    else
      value =	right_width/left_width;

    return value;
  }


	bool PeakShape::iteratorsSet() const
	{
		if(left_iterator_set_ && right_iterator_set_) return true;
		else return false;
	}

	PeakShape::PeakIterator PeakShape::getLeftEndpoint() const
	{
		return left_endpoint_;
	}
	void PeakShape::setLeftEndpoint(PeakShape::PeakIterator left_endpoint)
	{
		left_endpoint_=left_endpoint;
		left_iterator_set_ = true;
	}
	PeakShape::PeakIterator PeakShape::getRightEndpoint() const
	{
		return right_endpoint_;
	}

	void PeakShape::setRightEndpoint(PeakShape::PeakIterator right_endpoint)
	{
		right_endpoint_=right_endpoint;
		right_iterator_set_ = true;
	}

	
}
