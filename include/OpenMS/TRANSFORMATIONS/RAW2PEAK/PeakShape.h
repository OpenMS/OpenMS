// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKSHAPE_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKSHAPE_H

#include <cmath>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  /** 
  	@brief Internal representation of a peak shape (used by the PeakPickerCWT)

    It defines an asymmetric lorentzian and asymmetric hyperbolic squared secan function. 
  */
  struct OPENMS_DLLAPI PeakShape
  {
	  /** 
	    @brief Peak shape type (asymmetric lorentzian or asymmetric hyperbolic secans squared).
	
	    The peak shape can represent an asymmetric lorentzian function, given by 
	                  
	    l(x) = height/(1.+pow(left_width*(x - mz_position), 2)) (x<=mz_position) 
	                  
	    l(x) = height/(1.+pow(right_width*(x - mz_position), 2)) (x>mz_position)
	                  
	    or an asymmetric hyperbolic secans squared function 
	                  
	    s(x) = height/pow(cosh(left_width*(x-mz_position)), 2) (x<=mz_position)
	                  
	    s(x) = height/pow(cosh(right_width*(x-mz_position)), 2) (x>mz_position)
	  */
    enum Type
    {
      LORENTZ_PEAK,
      SECH_PEAK,
      UNDEFINED
    };
  
		/// Iterator to the raw data vector
		typedef MSSpectrum<>::const_iterator PeakIterator;

    /// Default constructor
    PeakShape()
        : height(0),
					mz_position(0),
					left_width(0),
					right_width(0),
					area(0),
					r_value(0),
					signal_to_noise(0.),
					type(UNDEFINED)
    {
			left_endpoint_ = exp_.end();
			right_endpoint_ = exp_.end();
    }
    
    /// Constructor that sets most of the members
    PeakShape(DoubleReal height_, DoubleReal mz_position_, DoubleReal left_width_, DoubleReal right_width_, DoubleReal area_, PeakIterator left_, PeakIterator right_, Type type_);

		/// Constructor that sets most of the members 
    PeakShape(DoubleReal height_, DoubleReal mz_position_, DoubleReal left_width_, DoubleReal right_width_, DoubleReal area_, Type type_);

		
    /// Copy constructor
    PeakShape(const PeakShape& rhs);
    
    /// Destructor
    virtual ~PeakShape()
    {
    }
    
    /// Assignment operator
    PeakShape& operator=(const PeakShape& rhs);
    
    //Equality operator
		bool operator==(const PeakShape& rhs) const;
    //Equality operator
		bool operator!=(const PeakShape& rhs) const;
		
    /// Compute the intensity of the peaks shape at position x
    DoubleReal operator() (DoubleReal x) const;
    /// Computes symmetry measure of the peak shape, which is corresponds to th ratio of the left and right width parameters.
    DoubleReal getSymmetricMeasure() const;
    /// Estimates the full width at half maximum.
    DoubleReal getFWHM() const;
		/// Check if endpoint iterators 
		bool iteratorsSet() const;
		
		PeakIterator getLeftEndpoint() const;
		void setLeftEndpoint(PeakIterator left_endpoint);
		
		PeakIterator getRightEndpoint() const;
		void setRightEndpoint(PeakIterator right_endpoint);
    /// Maximum intensity of the peak shape
    DoubleReal height;
    /// Centroid position
    DoubleReal mz_position;
    /// Left width parameter
    DoubleReal left_width;
    /// Right width parameter
    DoubleReal right_width;
    /// Area of the peak shape
    DoubleReal area;
    /** @brief Correlation coefficient.
      
      It represents the squared pearson correlation coefficient with the original data (0 <= r_value <= 1).
    */
    DoubleReal r_value;
    /// The signal to noise ratio at the mz_position
    DoubleReal signal_to_noise;
    
    ///peak shape type
    Type type;

    /**
    	 @brief  Comparison of mz_positions.
    */
    class OPENMS_DLLAPI PositionLess
    {
	    public:
	
      inline bool operator () (const PeakShape& a, const PeakShape& b)
      {
        return (a.mz_position < b.mz_position);
      }
      
    };
	protected:
		/// Left peak endpoint in the data
    PeakIterator left_endpoint_;
    /// Right peak endpoint in the data
    PeakIterator right_endpoint_;
		/// Needed for initialisation of endpoint iterators
		MSSpectrum<> exp_;
		/// flag if left endpoint iterator differs from default value
		bool left_iterator_set_;
		/// flag if left endpoint iterator differs from default value
		bool right_iterator_set_;
  };
} // namespace OpenMS

#endif
