	// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPMATCHERREGRESSION_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPMATCHERREGRESSION_H

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseMapMatcher.h>

#include <iostream>
#include <utility>

#include <gsl/gsl_fit.h>

#define DEBUG_MAPMATCHING
#undef DEBUG_MAPMATCHING

namespace OpenMS
{

  /**
		@brief Map matching using linear regression.
		
	*/
  template <typename ElementT = Feature >
  class MapMatcherRegression
		: public BaseMapMatcher<ElementT>
  {
	 public:

		/** @name Type definitions
		 */
		//@{

		typedef BaseMapMatcher<ElementT> Base;
		typedef typename Base::ElementType ElementType;
		typedef typename Base::ElementPairVector ElementPairVector;
		//@}

		/// Constructor
		MapMatcherRegression()
			: BaseMapMatcher<ElementType>()
		{}

		/// equality operator
		bool operator == (const MapMatcherRegression& rhs)
		{
			return (BaseMapMatcher<ElementType>::operator == (rhs) );
		}

		/// Destructor
		virtual ~MapMatcherRegression()
		{}

		/// estimates the transformation for each grid cell
		void estimateTransform()
		{
			ElementPairVector selection; // stores the pairs contained in the current cell

			for (typename ElementPairVector::iterator pair_iter = this->element_pairs_.begin();
					 pair_iter != this->element_pairs_.end();
					 ++pair_iter)
			{
				// check whether fulfills our quality requirement.
				if (pair_iter->getQuality() > this->min_quality_ )
				{
					selection.push_back(*pair_iter);
				}
			}

      // build arrays
			int num = selection.size();
			if (num > 2)
			{

				double* x = new double[num];
				double* y = new double[num];

				// loop over all dimensions //TODO
				for (UInt d=0; d<1;d++)
				{

					for (int i=0; i<num;i++)
					{
						x[i] = selection[i].getFirst().getPosition()[d];
						y[i] = selection[i].getSecond().getPosition()[d];
					}

					// estimate the transform for this dimension
					double slope, intercept, cov00, cov01, cov11, sumsq;

					gsl_fit_linear(x, 1, y, 1, num, &intercept, &slope, &cov00, &cov01, &cov11, &sumsq);

					// create the transform and save it in the cell
					this->grid_.setSlope(slope);
					this->grid_.setIntercept(intercept);
				} // end for (d)

				delete [] x;
				delete [] y;

			} // end for (gridcell)

		} // end void estimateTransform()

  }
  ; // end of class MapMatcherRegression

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_MAPMATCHERREGRESSION_H
