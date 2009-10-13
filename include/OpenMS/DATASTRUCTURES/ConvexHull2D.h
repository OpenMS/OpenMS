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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_CONVEXHULL2D_H
#define OPENMS_DATASTRUCTURES_CONVEXHULL2D_H

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

#include <vector>
#include <map>
#include <cmath>

#ifdef _MSC_VER // disable some CGAL warnings that distract from ours
#	pragma warning( push ) // save warning state
#	pragma warning( disable : 4244 )
#endif
#include <CGAL/Cartesian.h>
#include <CGAL/convex_hull_2.h>
#ifdef _MSC_VER
#	pragma warning( pop )  // restore old warning state
#endif


namespace OpenMS
{
	/**	
		@brief A 2-dimensional convex hull representation (conterclockwise)
		
		@ingroup Datastructures
	*/
	class OPENMS_DLLAPI ConvexHull2D
	{
		public:
			typedef DPosition<2> PointType;
			typedef std::vector< PointType > PointArrayType;
			typedef PointArrayType::size_type SizeType;
			typedef PointArrayType::const_iterator PointArrayTypeConstIterator;
			/// default constructor
			ConvexHull2D()
				: points_()
			{
			}

			/// Constructor from point array
			ConvexHull2D(const PointArrayType& points)
			{
				*this = points;
			}
			
			/// assignment operator
			ConvexHull2D& operator=(const ConvexHull2D& rhs)
			{
				if (&rhs == this) return *this;
				
				points_ = rhs.points_;
				
				return *this;
			}
			

			///assignment operator from a vector of points
			ConvexHull2D& operator=(const PointArrayType& points)
			{
				//init
				points_.clear();
				if (points.size()==0) return *this;
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
				return *this;
			}

			/// equality operator
			bool operator==(const ConvexHull2D& rhs) const
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
			DBoundingBox<2> getBoundingBox() const
			{
				DBoundingBox<2> bb;
				
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
					//convert input to cgal
					std::vector<Point_2> cgal_points;
					for (PointArrayTypeConstIterator it = points_.begin(); it!=points_.end(); ++it)
		      {
						cgal_points.push_back( Point_2((*it)[0], (*it)[1]) );    
		      }
					//calculate convex hull
					std::vector<Point_2> cgal_result;
		  		CGAL::convex_hull_2( cgal_points.begin(), cgal_points.end(), std::inserter(cgal_result, cgal_result.begin() ) );
		      // add the points
		      points_.clear();
					for (std::vector<Point_2>::const_iterator cit = cgal_result.begin(); cit !=	cgal_result.end(); ++cit)
					{
						points_.push_back( PointType(cit->x(),cit->y()) );			
					}
					return true;
				}
				return false;
			}
			
			/** 
				@brief returns if the @p point lies in the convex hull
			
			   @bug Due to numerical instabilities, this test for inclusion might return false even if
			        @p point is included  in the convex hull. As a preliminary workaround, we re-compute
			        the convex hull. In the meantime the points_ variable
			        is made mutable for the workaround. Undo that as soon as it is fixed. (Clemens)
			**/
			bool encloses(const PointType& point) const
			{
				if (!getBoundingBox().encloses(point))
				{
					return false;
				}				
				// re-calculate convex hull
				std::vector<Point_2> cgal_points;
				std::vector<Point_2> cgal_result;
				
				for (PointArrayTypeConstIterator it = points_.begin(); it!=points_.end(); ++it)
	      {
					cgal_points.push_back( Point_2((*it)[0], (*it)[1]) );    
	      }
				points_.clear();
					
				CGAL::convex_hull_2( cgal_points.begin(), cgal_points.end(), std::inserter(cgal_result, cgal_result.begin() ) );
	      // add the points
				for (std::vector<Point_2>::const_iterator cit = cgal_result.begin(); cit !=	cgal_result.end(); ++cit)
				{
					points_.push_back( PointType(cit->x(),cit->y()) );			
				} 					
				
				//convert input to cgal
				cgal_points.clear();
				cgal_result.clear();
				for (PointArrayTypeConstIterator it = points_.begin(); it!=points_.end(); ++it)
	      {
					cgal_points.push_back( Point_2((*it)[0], (*it)[1]) );    
	      }
				// add point to be tested
	      cgal_points.push_back( Point_2(point[0],point[1]) ); 
				
				//calculate convex hull
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
	      return true;
			}
			
		protected:
			mutable PointArrayType points_;
			typedef CGAL::Cartesian<DoubleReal>::Point_2 Point_2;
	};
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DCONVEXHULL_H
