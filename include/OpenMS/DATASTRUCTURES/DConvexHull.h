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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DCONVEXHULL_H
#define OPENMS_DATASTRUCTURES_DCONVEXHULL_H

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

#include <vector>
#include <map>
#include <math.h>

#ifdef CGAL_DEF
	#include <CGAL/Cartesian.h>
	#include <CGAL/convex_hull_2.h>
#endif

namespace OpenMS
{
	/**	
		@brief A d-dimensional convex hull representation (conterclockwise)
		
		@note If OpenMS is compiled with CGAL, the convex hull is calculated using CGAL. 
					Otherwise the gift-wrapping algorithm is used.
		
		@ingroup Datastructures
	*/
	template <Size D>
	class DConvexHull
	{
		public:
			typedef DPosition<D> PointType;
			typedef std::vector< DPosition<D> > PointArrayType;
			typedef typename PointArrayType::size_type SizeType;
			typedef typename PointArrayType::const_iterator PointArrayTypeConstIterator;
			/// default constructor
			DConvexHull()
				: points_()
			{
				if (D!=2)
				{
					throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
			}
			
			/// assignment operator
			DConvexHull& operator=(const DConvexHull& rhs)
			{
				if (&rhs == this) return *this;
				
				points_ = rhs.points_;
				
				return *this;
			}
			

			///constructor from a vector of points
			DConvexHull& operator=(const PointArrayType& points)
			{
				//init
				points_.clear();
				if (points.size()==0) return *this;
#ifdef CGAL_DEF				
				//convert input to cgal
				std::vector<Point_2> cgal_points;
				for (PointArrayTypeConstIterator it = points.begin(); it!=points.end(); ++it)
	      {
					cgal_points.push_back( Point_2((*it)[0], (*it)[1]) );    
	      }
				//calculate convex hull
				std::vector<Point_2> cgal_result;
	  		CGAL::convex_hull_2( cgal_points.begin(), cgal_points.end(), std::inserter(cgal_result, cgal_result.begin() ) );
	      // add the points
				for (std::vector<Point_2>::const_iterator cit = cgal_result.begin(); cit !=	cgal_result.end(); ++cit)
				{
					points_.push_back( PointType(cit->x(),cit->y()) );			
				} 	
#else
				giftWrapping_(points,points_);
#endif
				return *this;
			}

			/// equality operator
			bool operator==(const DConvexHull& rhs) const
			{
				// different size => return false
				if (points_.size()!= rhs.points_.size()) return false;
				
				//different points now => return false
				for (PointArrayTypeConstIterator it = rhs.points_.begin(); it !=	rhs.points_.end(); ++it)
				{
					if (find(points_.begin(),points_.end(),*it)==points_.end())
					{
						return false;
					}
				}
				return true;
			}
			
			/// removes all points
			void clear()
			{
				points_.clear();
			}

			/// accessor for the points
			const PointArrayType& getPoints() const
			{
				return points_;
			}
			
			/// returns the bounding box of the convex hull points
			DBoundingBox<D> getBoundingBox() const
			{
				DBoundingBox<D> bb;
				
				for (PointArrayTypeConstIterator it = points_.begin(); it!=points_.end(); ++it)
				{
					bb.enlarge(*it);
				}
				
				return bb;
			}
			
			/// adds a point to the convex hull if it is not already contained. Returns if the point was added.
			bool addPoint(const PointType& point)
			{
				if (!encloses(point))
				{
					points_.push_back(point);
					return true;
				}
				return false;
			}
			
			/// returns if the @p point lies in the convex hull
			bool encloses(const PointType& point) const
			{
				if (!getBoundingBox().encloses(point))
				{
					return false;
				}
#ifdef CGAL_DEF				
				//convert input to cgal
				std::vector<Point_2> cgal_points;
				for (PointArrayTypeConstIterator it = points_.begin(); it!=points_.end(); ++it)
	      {
					cgal_points.push_back( Point_2((*it)[0], (*it)[1]) );    
	      }
	      cgal_points.push_back( Point_2(point[0],point[1]) ); 
				
				//calculate convex hull
				std::vector<Point_2> cgal_result;
	  		CGAL::convex_hull_2( cgal_points.begin(), cgal_points.end(), std::inserter(cgal_result, cgal_result.begin() ) );
				
				//point added => return false
				if (cgal_result.size()!=points_.size()) return false;
				
				//different points now => return false
				for (std::vector<Point_2>::const_iterator cit = cgal_result.begin(); cit !=	cgal_result.end(); ++cit)
				{
					if (find(points_.begin(),points_.end(),PointType(cit->x(),cit->y()))==points_.end())
					{
						return false;
					}
				}
#else
			PointArrayType tmp = points_,new_hull;
			tmp.push_back(point);
			giftWrapping_(tmp,new_hull);
			
			//point added => return false
			if (new_hull.size()!=points_.size()) return false;
			
			//different points now => return false
			for (PointArrayTypeConstIterator it = new_hull.begin(); it !=	new_hull.end(); ++it)
			{
				if (find(points_.begin(),points_.end(),*it)==points_.end())
				{
					return false;
				}
			}
#endif
	      return true;
			}
			
		protected:
			PointArrayType points_;
#ifdef CGAL_DEF
			typedef CGAL::Cartesian<double>::Point_2 Point_2;
#endif
		
		/// fallback implementation if CGAL is not present
		void giftWrapping_(const PointArrayType& input, PointArrayType& output) const
		{
			//nothing to do for one or two points
			if (input.size()<3)
			{
				output = input;
				return;
			}
			else
			{
				output.clear();
			}
			//Precision for double comparisons
			const double PRECISION = 0.0001;
			
			// keep track of already in hull included peaks to avoid unnecessary computations of triangle area
			std::map<PointType, bool> is_included;
			
			typename PointType::CoordinateType min_mz = std::numeric_limits<typename PointType::CoordinateType>::max();
			PointArrayTypeConstIterator min = input.begin();
			
			// Find peak with minimal mz to start wrapping
			for (PointArrayTypeConstIterator it = input.begin(); it!=input.end(); ++it)
			{
				if ((*it)[1] < min_mz)
				{
				   min_mz = (*it)[1];
				   min = it;
				}
				is_included[*it] = false;
			}
			output.push_back(*min);
			
			// Hull peaks denoting current hull line
			PointArrayTypeConstIterator hull_peak1 = min;
			PointArrayTypeConstIterator start = input.begin();
			if (start==min)
			++start;  // don't start at "min" because of while-condition
			PointArrayTypeConstIterator hull_peak2 = start;
			
			while (hull_peak2!=min)
			{
				bool found_any = false;
				for (PointArrayTypeConstIterator it = input.begin(); it!=input.end(); ++it)
				{
					// skip if already used
					if (is_included[*it] || it==hull_peak1 || it==hull_peak2)
					   continue;
					
					found_any = true;
					// "it" lies to the right of the line [hull_peak1,hull_peak2]
					
					// triangle area via determinant: x0*y1+x1*y2+x2*y0-x2*y1-x1*y0-x0*y2
					double area = (*hull_peak1)[1]*(*hull_peak2)[0] + (*hull_peak2)[1]*(*it)[0] + (*it)[1]*(*hull_peak1)[0] - (*it)[1]*(*hull_peak2)[0] - (*hull_peak2)[1]*(*hull_peak1)[0] - (*hull_peak1)[1]*(*it)[0];

					if (area>-PRECISION)
					{
						// area almost 0 -> collinear input
						// -> avoid consecutive peaks with equal mz or rt coordinate
						if (fabs(area)<PRECISION)
						{
						   double mz1 = (*hull_peak1)[1];
						   double mz2 = (*hull_peak2)[1];
						   double mz3 = (*it)[1];
						   double rt1 = (*hull_peak1)[0];
						   double rt2 = (*hull_peak2)[0];
						   double rt3 = (*it)[0];
						   if ( 
						   			(fabs(mz2-mz3)<PRECISION && fabs(rt2-rt1) > fabs(rt3-rt1)) 
						   			|| 
						   			(fabs(rt2-rt3)<PRECISION && fabs(mz2-mz1) > fabs(mz3-mz1))
						   			)
						   {
						       is_included[*it] = true;
						       continue;
						   }
						}
						hull_peak2 = it;  // "it" becomes new hull peak
			  	}
				}
				
				if (!found_any)
				{
				   hull_peak2 = min; // no available peaks anymore
				   continue;
				}
				
				if (hull_peak2 == min)
				{
					continue;  // finish loop
				}
				is_included[*hull_peak2] = true;
				
				// continue wrapping
				hull_peak1 = hull_peak2;
				// hull_peak2 satisfies the contition: all peaks lie to the left of [hull_peak1,hull_peak2]
				output.push_back(*hull_peak2);
				
				start = input.begin();
				if (start==min)
				{
				   ++start;  // don't start at "min" because of while-condition
				}
				hull_peak2 = start;
			}
		}
	};
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DCONVEXHULL_H
