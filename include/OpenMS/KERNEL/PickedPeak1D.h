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

#ifndef OPENMS_KERNEL_PICKEDPEAK1D_H
#define OPENMS_KERNEL_PICKEDPEAK1D_H

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShapeType.h>

#include <math.h>

namespace OpenMS
{
  /**
  	@brief 1-dimensional picked peak.
  	
  	Extends Peak1D by members for peak-picking algorithms.
  	
  	The intensity of a peak is defined as the maximum of the model fitted to the raw data during peak picking
  	i.e. aproximately the height of the highest raw data point.
  	
  	@ingroup Kernel
  */
  class PickedPeak1D
    : public Peak1D
  {
  public:

    /// Default constructor
    PickedPeak1D():
        Peak1D(),
        r_value_(0),
        area_(0),
        fwhm_(0),
        left_width_paramter_(0),
        right_width_paramter_(0),
        type_(PeakShapeType::UNDEFINED),
        charge_(0),
        signal_to_noise_(0)
    {
    }

    /// Copy sontructor
    inline PickedPeak1D(PickedPeak1D const& p)
      :
        Peak1D(p),
        r_value_(p.r_value_),
        area_(p.area_),
        fwhm_(p.fwhm_),
        left_width_paramter_(p.left_width_paramter_),
        right_width_paramter_(p.right_width_paramter_),
        type_(p.type_),
        charge_(p.charge_),
        signal_to_noise_(p.signal_to_noise_)
    {
    }

    /// Desctructor
    ~PickedPeak1D()
    {
    }
    
    /// Non-mutable access to the correlation coefficient between raw data and the peak model
    inline const DoubleReal& getRValue() const { return r_value_; }
    /// Mutable access to the peak correlation coefficient between raw data and the peak model
    inline DoubleReal& getRValue() { return r_value_; }
    /// Mutable access to the peak correlation coefficient between raw data and the peak model
    inline void setRValue(const DoubleReal& r_value) { r_value_ = r_value; }

    /// Non-mutable access to the peak area
    inline const DoubleReal& getArea() const { return area_; }
    /// Mutable access to the peak area
    inline DoubleReal& getArea() { return area_; }
    /// Mutable access to the peak area
    inline void setArea(const DoubleReal& area) { area_ = area; }

    /// Non-mutable access to the peak FWHM
    inline const DoubleReal& getFWHM() const { return fwhm_; }
    /// Mutable access to the peak FWHM
    inline DoubleReal& getFWHM() { return fwhm_; }
    /// Mutable access to the peak FWHM
    inline void setFWHM(const DoubleReal& fwhm) { fwhm_ = fwhm; }

    /// Non-mutable access to the width parameter of the left peak side
    inline const DoubleReal& getLeftWidthParameter() const { return left_width_paramter_; }
    /// Mutable access to the width parameter of the left peak side
    inline DoubleReal& getLeftWidthParameter() { return left_width_paramter_; }
    /// Mutable access to the width parameter of the left peak side
    inline void setLeftWidthParameter(const DoubleReal& left_width_paramter) { left_width_paramter_ = left_width_paramter; }

    /// Non-mutable access to the width parameter of the right peak side
    inline const DoubleReal& getRightWidthParameter() const { return right_width_paramter_; }
    /// Mutable access to the width parameter of the right peak side
    inline DoubleReal& getRightWidthParameter() { return right_width_paramter_; }
    /// Mutable access to the width parameter of the right peak side
    inline void setRightWidthParameter(const DoubleReal& right_width_paramter) { right_width_paramter_ = right_width_paramter; }

    /// Non-mutable access to the peak shape
    inline const PeakShapeType::Enum& getPeakShape() const { return type_; }
    /// Mutable access to the peak shape
    inline PeakShapeType::Enum& getPeakShape() { return type_; }
    /// Mutable access to the peak shape
    inline void setPeakShape(const PeakShapeType::Enum& type) { type_ = type; }

    /// Non-mutable access to the peak charge
    inline const SignedInt& getCharge() const { return charge_; }
    /// Mutable access to the peak charge (Set to 0 if unknown)
    inline SignedInt& getCharge() { return charge_; }
    /// Mutable access to the peak charge (Set to 0 if unknown)
    inline void setCharge(const SignedInt& charge) { charge_ = charge; }

    /// Non-mutable access to the signal to noise value
    inline const DoubleReal& getSN() const { return signal_to_noise_; }
    /// Mutable access to the the signal to noise value
    inline DoubleReal& getSN() { return signal_to_noise_; }
    /// Mutable access to the the signal to noise value
    inline void setSN(const DoubleReal& signal_to_noise) { signal_to_noise_ = signal_to_noise; }
    
    /// Assignment operator
    PickedPeak1D& operator = (const PickedPeak1D& rhs)
    {
      if (this==&rhs) return *this;

      Peak1D::operator = (rhs);
      r_value_             = rhs.r_value_;
      area_                = rhs.area_;
      fwhm_                = rhs.fwhm_;
      left_width_paramter_ = rhs.left_width_paramter_;
      right_width_paramter_= rhs.right_width_paramter_;
      type_                = rhs.type_;
      charge_              = rhs.charge_;
      signal_to_noise_     = rhs.signal_to_noise_;

      return *this;
    }

    /// Equality operator
    bool operator == (const PickedPeak1D& rhs) const
    {
    	return r_value_ == rhs.r_value_ && 
    				 area_ == rhs.area_ &&	
    				 fwhm_ == rhs.fwhm_ &&  
    				 type_ == rhs.type_ && 
    				 charge_ == rhs.charge_  && 
    				 left_width_paramter_ == rhs.left_width_paramter_ && 
    				 right_width_paramter_ == rhs.right_width_paramter_ && 
    				 signal_to_noise_ == rhs.signal_to_noise_ && 
    				 Peak1D::operator==(rhs);
    }

    /// Equality operator
    bool operator != (const PickedPeak1D& rhs) const
    {
      return !(operator == (rhs));
    }

    /**
    	 @brief Returns the symmetry s of a peak with: (asymmetric peaks) 0 < s <= 1 (symmetric peaks).
    */
    double getSymmetricMeasure() const
    {
      double value=0.;

      if (left_width_paramter_ < right_width_paramter_)
        value = left_width_paramter_/right_width_paramter_;
      else
        value =	right_width_paramter_/left_width_paramter_;

      return value;
    }

    /**
    	 @brief Returns the value of the peak shape function at position x.
    */
    double operator () (const double& x) const
    {
      double value;

      switch (type_)
      {
      case PeakShapeType::LORENTZ_PEAK:
        if (x <= getMZ())
          value = getIntensity()/(1.+pow(left_width_paramter_*(x - getMZ()), 2));
        else
          value = getIntensity()/(1.+pow(right_width_paramter_*(x - getMZ()), 2));
        break;

      case PeakShapeType::SECH_PEAK:
        if ( x <= getMZ())
          value = getIntensity()/pow(cosh(left_width_paramter_*(x-getMZ())), 2);
        else
          value = getIntensity()/pow(cosh(right_width_paramter_*(x-getMZ())), 2);
        break;

      default:
        value = -1.;
        break;
      }

      return value;
    }

    /**	@name	Comparator classes.
    	These classes implement binary predicates that can be used 
    	to compare two peaks with respect to their width.
    	They are usually employed by the sort methods of PickedPeak1DList and PickedPeak1DArray.
    */
    //@{
    /**
    	 @brief Comparator for the width.

    	 Lexicographical comparison from dimension 0 to dimension D is done.
    */
    struct WidthLess
          : public std::binary_function <PickedPeak1D, PickedPeak1D, bool>
    {
      inline bool operator () (const PickedPeak1D& a, const PickedPeak1D& b)
      {
        double a_width=0;
        double b_width=0;

        switch (a.type_)
        {
        case PeakShapeType::LORENTZ_PEAK:
          {
            a_width = sqrt(10*a.getIntensity()-1)/a.left_width_paramter_;
            a_width+= sqrt(10*a.getIntensity()-1)/a.right_width_paramter_;
          }
          break;

        case PeakShapeType::SECH_PEAK:
          {
            a_width = acosh(sqrt(3*a.getIntensity())/3)/a.left_width_paramter_;
            a_width+= acosh(sqrt(3*a.getIntensity())/3)/a.right_width_paramter_;
          }
          break;

        default:
          {
            a_width = -1.;
          }
          break;
        }

        switch (b.type_)
        {
        case PeakShapeType::LORENTZ_PEAK:
          {
            b_width = sqrt(10*b.getIntensity()-1)/b.left_width_paramter_;
            b_width+= sqrt(10*b.getIntensity()-1)/b.right_width_paramter_;
          }
          break;

        case PeakShapeType::SECH_PEAK:
          {
            b_width = acosh(sqrt(3*b.getIntensity())/3)/b.left_width_paramter_;
            b_width+= acosh(sqrt(3*b.getIntensity())/3)/b.right_width_paramter_;
          }
          break;

        default:
          {
            b_width= -1.;
          }
          break;
        }

        return (a_width < b_width);
      }

      /**
      	 @brief Operator to check if comparison is done increasing or decreasing.
      	
      	 Sometimes we need a way to find out which way the CoordinateType is
      	 sorted and adding this overload seems to be the best way to achieve that goal.
      */
      inline bool operator () ( PositionType const & left, PositionType const & right ) const throw()
      {
        return (left < right );
      }
    };
    //@}

  protected:
    /// The correlation factor (degree how the raw data peak matches with a computed peak (lorentzian or sech)
    DoubleReal r_value_;
    /// The area
    DoubleReal area_;
    /// Full-width-at-half-max
    DoubleReal fwhm_;
    /// The function dependent left and right width parameter
    DoubleReal	left_width_paramter_;
    DoubleReal	right_width_paramter_;
    /// The function that was used for fitting the peak shape
    PeakShapeType::Enum type_;
    /// The peak charge
    SignedInt charge_;
    /// The signal to noise value of the peak
    DoubleReal signal_to_noise_;

  };

} // namespace OpenMS

#endif // OPENMS_KERNEL_DPICKEDPEAK_H
