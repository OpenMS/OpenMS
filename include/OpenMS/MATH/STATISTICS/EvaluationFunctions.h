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
// $Id: EvaluationFunctions.h,v 1.7 2006/04/12 08:21:18 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_STATISTICS_EVALUATIONFUNCTIONS_H
#define OPENMS_MATH_STATISTICS_EVALUATIONFUNCTIONS_H

#include<OpenMS/CONCEPT/Types.h>

#include <cmath>

namespace OpenMS 
{
  /**
    @brief Serves for calculation of result statistics
    
    Offers functions like meanSquareError or pearsonCorrelationCoefficient

  	@todo add to Math namespace! merge with BasicStatistics? (Nico)
    
    @ingroup Math
  */
  class EvaluationFunctions
  {
    public:
      
      template <class iteratorT>
      static DoubleReal meanSquareError(iteratorT begin_a, 
      																 iteratorT end_a, 
      																 iteratorT begin_b, 
      																 iteratorT end_b)
      {
      	SignedInt count = 0;
      	DoubleReal error = 0;
      	iteratorT it_a = begin_a;
      	iteratorT it_b = begin_b;
      	      	
      	while(it_a != end_a && it_b != end_b)
      	{
      		error += pow(*it_a - *it_b, 2);
      		count++;
      		it_a++;
      		it_b++;      		
      	}
      	
      	return error / count;
      	
      }

      template <class iteratorT>
      static DoubleReal pearsonCorrelationCoefficient(iteratorT begin_a, 
      													 										 iteratorT end_a, 
      													 										 iteratorT begin_b, 
      													 										 iteratorT end_b)
      {
      	SignedInt count = 0;
      	DoubleReal sum_a = 0;
      	DoubleReal sum_b = 0;
      	DoubleReal mean_a = 0;
      	DoubleReal mean_b = 0;      	
      	iteratorT it_a = begin_a;
      	iteratorT it_b = begin_b;      	
      	DoubleReal numerator = 0;
      	DoubleReal denominator_a = 0;
      	DoubleReal denominator_b = 0;
      	DoubleReal temp_a;
      	DoubleReal temp_b;
      	      	
      	while(it_a != end_a && it_b != end_b)
      	{
					sum_a += *it_a;
					sum_b += *it_b;
      		count++;
      		it_a++;
      		it_b++;      		
      	}
      	mean_a = sum_a / count;
      	mean_b = sum_b / count;
      	
      	it_a = begin_a;
      	it_b = begin_b;
      	while(it_a != end_a && it_b != end_b)
      	{
      		temp_a = *it_a - mean_a;
      		temp_b = *it_b - mean_b;
      		
      		numerator += (temp_a * temp_b);
      		denominator_a += pow(temp_a, 2);
      		denominator_b += pow(temp_b, 2);
      		it_a++;
      		it_b++;      		
				}      	
      	
      	return numerator / sqrt(denominator_a * denominator_b);      	
      }
    
    private:

      /// Constructor
      EvaluationFunctions();
      /// Copy constructor
      EvaluationFunctions(const EvaluationFunctions& source);
      /// Destructor
      ~EvaluationFunctions();
      
      /// Assignment operator
      EvaluationFunctions& operator = (const EvaluationFunctions& source);

  };
 
} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_EVALUATIONFUNCTIONS_H
