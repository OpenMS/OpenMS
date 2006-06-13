// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: ChiSquared.h,v 1.4 2006/03/28 16:19:59 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_STATISTICS_CHISQUARED_H
#define OPENMS_MATH_STATISTICS_CHISQUARED_H

#include <vector>

#include <gsl/gsl_cdf.h>

namespace OpenMS
{
  /**
  	@brief Chi squared computation, includes binning.

  	@todo add to Math namespace (Clemens)
  	
  	@ingroup Math
  */
  template < typename Coordinate = double, typename Probability = Coordinate >
  class
  ChiSquared
  {

	 public:

    typedef Probability probability_type;
    typedef Coordinate coordinate_type;
    typedef std::vector < probability_type > container_type;

    /// Default constructor.  You will need to call calculateChiSquared() manually.
    ChiSquared () { reset(); }

    // No particular copy constructor needed (?)

    /// This does the actual calculation.
    void
    calculateChiSquared ()
    {
      int rightMost = model_->size() - 1;
      int lastLeft = 0;
      int currentLeft = 0;
      int lastRight = rightMost;
      int currentRight = rightMost;

      chiSquared_ = 0;

      probability_type quantileObservation;
      probability_type quantileModel;
      probability_type diff;

#define VERBOSEBINNING 0
      for (;;) {
				//---
				quantileModel = 0;
				quantileObservation = 0;
				do {
					quantileModel += (*model_)[currentLeft];
					quantileObservation += (*observation_)[currentLeft];
#if VERBOSEBINNING
					DUMP(currentLeft);
#endif
					++currentLeft;
				} while ( quantileModel < threshold_ && currentLeft < currentRight );
				diff = quantileModel - quantileObservation;
				diff *= diff;
				chiSquared_ += diff / quantileModel;
				++numberOfBins_;
#if VERBOSEBINNING
				cerr << lastLeft << ' ' << currentLeft << ' ' << currentRight << ' ' << lastRight << ' ' << rightMost << ' ' 
						 << quantileModel << ' ' << quantileObservation << ' ' << diff << ' ' << chiSquared_ << ' ' << numberOfBins_ << "\n";
#endif
				if ( currentLeft > currentRight ) break;
				lastLeft = currentLeft;
				//---
				quantileModel = 0;
				quantileObservation = 0;
				do {
					quantileModel += (*model_)[currentRight];
					quantileObservation += (*observation_)[currentRight];
#if VERBOSEBINNING
					DUMP(currentRight);
#endif
					--currentRight;
				} while ( quantileModel < threshold_ && currentLeft < currentRight );
				diff = quantileModel - quantileObservation;
				diff *= diff;
				chiSquared_ += diff / quantileModel;
				++numberOfBins_;
#if VERBOSEBINNING
				cerr << lastLeft << ' ' << currentLeft << ' ' << currentRight << ' ' << lastRight << ' ' << rightMost << ' ' 
						 << quantileModel << ' ' << quantileObservation << ' ' << diff << ' ' << chiSquared_ << ' ' << numberOfBins_ << "\n";
#endif
				if ( currentLeft > currentRight ) break;
				lastRight = currentRight;
				//---
      }
      return;
    }
#undef VERBOSEBINNING


    void
    reset ()
    {
      model_ = 0;
      observation_ = 0;
      threshold_ = 0;
      chiSquared_ = 0;
      numberOfBins_ = 0;
    }


    void
    setModel ( container_type const & _model )
    { model_ = &_model; }

    container_type &
    getModel ()
    { return *model_; }


    void
    setObservation ( container_type const & _observation )
    { observation_ = &_observation; }

    container_type &
    getObservation ()
    { return *observation_; }


    void
    setThreshold ( probability_type const _threshold )
    { threshold_=_threshold; return; }

    probability_type
    getThreshold()
    { return threshold_; }


    probability_type
    getChiSquared ()
    { return chiSquared_; }


    int
    getNumberOfBins ()
    { return numberOfBins_; }


    double
    getGoodnessOfFit( int degreesOfFreedom )
    {   
      return getGoodnessOfFit ( chiSquared_, degreesOfFreedom );
    }

    static double
    getGoodnessOfFit( probability_type _chiSquared, int degreesOfFreedom )
    {   
      return gsl_cdf_chisq_Q ( _chiSquared, degreesOfFreedom );
    }


	 protected:

    container_type const * model_;
    container_type const * observation_;

    probability_type threshold_;
    probability_type chiSquared_;

    int numberOfBins_;

  };
  
} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_CHISQUARED_H
