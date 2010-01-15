// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_STATISTICS_ROCCURVE_H
#define OPENMS_MATH_STATISTICS_ROCCURVE_H

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>

#ifdef _MSC_VER // disable some CGAL warnings that distract from ours
#	pragma warning( push ) // save warning state
#	pragma warning( disable : 4244 )
#endif
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#ifdef _MSC_VER
#	pragma warning( pop )  // restore old warning state
#endif

#include <list>
#include <vector>

namespace OpenMS
{
	namespace Math
	{
	  /**
	  	@brief ROCCurves show the tradeoff in sensitivity and specitivity for binary classifiers using different cutoff values
	  	
	  	@ingroup Math
	  */
	  class OPENMS_DLLAPI ROCCurve
	  {
	  public:
		
	    // @name Constructors and Destructors
	    // @{
			/// default constructor
	    ROCCurve();
	
			/// destructor
	    virtual ~ROCCurve();
	
			/// copy constructor
	    ROCCurve(const ROCCurve& source);
			// @}
	
			// @name Operators
			// @{
			/// assignment operator
	    ROCCurve& operator = (const ROCCurve& source);
	    // @}
	
			// @name Accessors
			// @{
	    /// insert score, type pair 
	    void insertPair(double score, bool clas);
	
	    /// returns Area Under Curve
	    double AUC();
	
	    /// some points in the ROC Curve
	    std::vector<std::pair<double, double> > curve(UInt resolution = 10);
	
			///
	    double cutoffPos(double fraction = 0.95);
	
			///
	    double cutoffNeg(double fraction = 0.95);
			// @}
	
	  private:

      typedef CGAL::Point_2<CGAL::Cartesian<double> > Point;
      typedef CGAL::Polygon_2<CGAL::Cartesian<double> > Polygon;

      /// predicate for sort()
      class OPENMS_DLLAPI simsortdec
      {
      	public:
        	
					bool operator () (const std::pair<double,bool>& a, const std::pair<double,bool>& b)
        	{
          	return b.first < a.first;
        	}
      };
			

	    std::list<std::pair<double,bool> > score_clas_pairs_;
			
	    UInt pos_;
			
	    UInt neg_;
	  };
	}
}
#endif // OPENMS_MATH_STATISTICS_ROCCURVE_H
