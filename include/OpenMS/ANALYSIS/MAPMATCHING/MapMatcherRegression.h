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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DMAPMATCHERREGRESSION_H
#define OPENMS_ANALYSIS_MAPMATCHING_DMAPMATCHERREGRESSION_H

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

		/// Copy constructor
		MapMatcherRegression(const MapMatcherRegression& source)
			: BaseMapMatcher<ElementType>(source)
		{}

		///  Assignment operator
		MapMatcherRegression& operator = (const MapMatcherRegression& source)
		{
			if (&source==this)
				return *this;

			BaseMapMatcher<ElementType>::operator = (source);
			return *this;
		}

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
			for ( typename Grid::iterator grid_iter = this->grid_.begin();
						grid_iter != this->grid_.end();
						++grid_iter
					)
			{
				ElementPairVector selection; // stores the pairs contained in the current cell

#ifdef DEBUG_MAPMATCHING

				std::cout << "Estimate the transformation with : " << this->element_pairs_.size() << " pairs " << std::endl;
#endif

				for (typename ElementPairVector::iterator pair_iter = this->element_pairs_.begin();
						 pair_iter != this->element_pairs_.end();
						 ++pair_iter)
				{
					// check whether the current feature is contained in the cell
					// and fulfills our quality requirement.
					if (grid_iter->encloses(pair_iter->getFirst().getPosition())
							&& pair_iter->getQuality() > this->min_quality_ )
					{
						selection.push_back(*pair_iter);

#ifdef DEBUG_MAPMATCHING
// 						std::cout << "Pair " << pair_iter->first.getPosition() << " " << pair_iter->second.getPosition() << std::endl;
#endif

					}
				} // end for (pair_iter)

	  // build arrays
				int num = selection.size();
				if (num > 2)
				{

					double* x = new double[num];
					double* y = new double[num];

					grid_iter->getMappings().clear();

					// loop over all dimensions
					for (UInt d=0; d<2;d++)
					{

						for (int i=0; i<num;i++)
						{
							x[i] = selection[i].getFirst().getPosition()[d];
							y[i] = selection[i].getSecond().getPosition()[d];
						}

						// estimate the transform for this dimension
						double slope, intercept, cov00, cov01, cov11, sumsq;

						gsl_fit_linear(x, 1, y, 1, num, &intercept, &slope, &cov00, &cov01, &cov11, &sumsq);

#ifdef DEBUG_MAPMATCHING

						std::cout << "Estimating transform for dimension " << d << std::endl;
						std::cout << "Best fit: Y = " << intercept << " + " << slope << "* X" << std::endl;
            std::cout << "Sumsquares: " << sumsq << std::endl;
            std::cout << "Pearson: " << (cov01*cov01)/(cov00*cov11) << std::endl;
#endif

						// create the transform and save it in the cell
						grid_iter->getMappings().push_back(new LinearMapping(slope,intercept));

					} // end for (d)

					delete [] x;
					delete [] y;

				} // end for (gridcell)
			}

		} // end void estimateTransform()

  }
  ; // end of class MapMatcherRegression

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DMAPMATCHERREGRESSION_H
