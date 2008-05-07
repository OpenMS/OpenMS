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
//  version 2.1 of the License, or (at your option) any later version
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMPOSECLUSTERING_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMPOSECLUSTERING_H

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>

#include <gsl/gsl_fit.h>

namespace OpenMS
{
	/**
		@brief A map alignment algorithm based on spectrum similarity (dynamic programming). 		
	*/
	class MapAlignmentAlgorithmPoseClustering
	 : public MapAlignmentAlgorithm
	{
		public:
			/// Default constructor
			MapAlignmentAlgorithmPoseClustering();

			/// Destructor
			virtual ~MapAlignmentAlgorithmPoseClustering();

			//Docu in base class
			virtual void alignPeakMaps(std::vector< MSExperiment<> >&);
				
			//Docu in base class
			virtual void alignFeatureMaps(std::vector< FeatureMap<> >&);
			
			///Creates a new instance of this class (for Factory)
			static MapAlignmentAlgorithm* create()
			{
				return new MapAlignmentAlgorithmPoseClustering();
			}
			
			///Returns the product name (for the Factory)
			static String getProductName()
			{
				return "pose_clustering";
			}
			
		private:

			///Copy constructor is not implemented -> private
			MapAlignmentAlgorithmPoseClustering(const MapAlignmentAlgorithmPoseClustering& );
			///Assignment operator is not implemented -> private
			MapAlignmentAlgorithmPoseClustering& operator=(const MapAlignmentAlgorithmPoseClustering& );
			
			template<typename ElementPairVector>
			LinearMapping calculateRegression_(const ElementPairVector& pairs)
			{
				UInt size = pairs.size();
				
				//create datastructures for GSL linear fit
				double* x = new double[size];
				double* y = new double[size];

				for (UInt i=0; i<size;i++)
				{
					x[i] = pairs[i].getFirst().getPosition()[0];
					y[i] = pairs[i].getSecond().getPosition()[0];
				}

				// estimate the transformation
				double slope, intercept, cov00, cov01, cov11, sumsq;
				gsl_fit_linear(x, 1, y, 1, size, &intercept, &slope, &cov00, &cov01, &cov11, &sumsq);
				
				//release memory
				delete [] x;
				delete [] y;

				//return result
				LinearMapping lm;
				lm.setSlope(slope);
				lm.setIntercept(intercept);
				return lm;
			}

	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMPOSECLUSTERING_H
