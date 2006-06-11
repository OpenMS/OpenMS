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
// $Id: DMapMatcherRegression.h,v 1.12 2006/04/18 15:22:54 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DMAPMATCHERREGRESSION_H
#define OPENMS_ANALYSIS_MAPMATCHING_DMAPMATCHERREGRESSION_H

#include<OpenMS/KERNEL/DFeature.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DLinearMapping.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DBaseMapMatcher.h>

#include<iostream.h>  
#include <utility>

#include<gsl/gsl_fit.h>

namespace OpenMS
{
	
	/**
	 	@brief Map matching using linear regression.
	 **/
	template <UnsignedInt D>		
	class DMapMatcherRegression
		: public DBaseMapMatcher<D>
	{
		public:
		
		/** @name Type definitions
		*/
		//@{	
		/// The grid is simply a vector of cells.
		typedef std::vector<DGridCell<2> > Grid;
		/// The feature pairs are computed by the feature matching class
		typedef std::vector<DFeaturePair<2> > FeaturePairVector;
		//@}
				
		/// Constructor
		DMapMatcherRegression()
			: DBaseMapMatcher<D>()  {}
		
		/// Copy constructor
		DMapMatcherRegression(const DMapMatcherRegression& source)
			: DBaseMapMatcher<D>(source) {}
				
		///  Assignment operator
    DMapMatcherRegression& operator = (const DMapMatcherRegression& source)
    {
    	if (&source==this) return *this;
    	
			DBaseMapMatcher<D>::operator = (source);
			return *this;
		}
		
		/// equality operator
		bool operator == (const DMapMatcherRegression& rhs)
		{
			return (DBaseMapMatcher<D>::operator == (rhs) );
		}
				
		/// Destructor
		virtual ~DMapMatcherRegression() {}
				
    /// estimates the transformation for each grid cell
    void estimateTransform()
		{
			for ( Grid::iterator grid_iter = this->grid_.begin();
						grid_iter != this->grid_.end();
						++grid_iter
					)
			{
								
				FeaturePairVector selection; // stores the pairs contained in the current cell
			
				for (FeaturePairVector::iterator pair_iter = this->feature_pairs_.begin();
						 pair_iter != this->feature_pairs_.end();
						 ++pair_iter)
				{
							
					// check whether the current feature is contained in the cell
					// and fulfills our quality requirement.
					if (grid_iter->encloses(pair_iter->getFirst().getPosition())
							&& pair_iter->getQuality() > this->min_quality_ )
					{
							selection.push_back(*pair_iter);
					}
					
				} // end for (pair_iter)
			
				// build arrays
				int num = selection.size();
			
				if (num == 0) 	continue;
								
				double* x = new double[num];
				double* y = new double[num];
			
				// loop over all dimensions 
				for (UnsignedInt d=0; d<D;d++)
				{
					
					for (int i=0; i<num;i++)
					{
						x[i] = selection[i].getFirst().getPosition()[d];
						y[i] = selection[i].getSecond().getPosition()[d];
						#ifdef DEBUG_MAPMATCHING
						std::cout << "x[" << i << "] = " << x[i] << std::endl;
						std::cout << "y[" << i << "] = " << y[i] << std::endl;
						#endif
					} 
				
					// estimate the transform for this dimension
					double slope, intercept, cov00, cov01, cov11, sumsq;
				
					gsl_fit_linear(x, 1, y, 1, num, &intercept, &slope, &cov00, &cov01, &cov11, &sumsq);        		

					#ifdef DEBUG_MAPMATCHING
					std::cout << "Estimating transform for dimension " << d << std::endl;
					std::cout << "Best fit: Y = " << intercept << " + " << slope << "* X" << std::endl; 
					std::cout << "Sumsquares: " << sumsq << std::endl;
					#endif
										
					// create the transform and save it in the cell
					grid_iter->getMappings().push_back(new DLinearMapping<1>(slope,intercept));
										
				} // end for (d)
				
				delete [] x;
				delete [] y;

			} // end for (gridcell)
			
		} // end void estimateTransform()   
								
	}; // end of class MapMatcherRegression
	
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DMAPMATCHERREGRESSION_H
