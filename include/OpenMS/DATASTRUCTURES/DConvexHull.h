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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DCONVEXHULL_H
#define OPENMS_DATASTRUCTURES_DCONVEXHULL_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

#include <vector>

#include <CGAL/Cartesian.h>
#include <CGAL/convex_hull_2.h>

namespace OpenMS
{
	/**	
		@brief A d-dimensional convex hull representation (conterclockwise)
		
		@ingroup Datastructures
	*/
	template <Size D>
	class DConvexHull
	{
		public:
			typedef DPosition<D> PointType;
			typedef std::vector< DPosition<D> > PointArrayType;
			typedef typename PointArrayType::size_type SizeType;
				
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
			

			/// constructor from a vector of points
			DConvexHull& operator=(const PointArrayType& points)
			{
				//init
				points_.clear();
				if (points.size()==0) return *this;
				
				//convert input to cgal
				std::vector<Point_2> cgal_points;
				for (typename PointArrayType::const_iterator it = points.begin(); it!=points.end(); ++it)
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
			bool operator==(const DConvexHull& rhs) const
			{
				// different size => return false
				if (points_.size()!= rhs.points_.size()) return false;
				
				//different points now => return false
				for (typename PointArrayType::const_iterator it = rhs.points_.begin(); it !=	rhs.points_.end(); ++it)
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
				
				for (typename PointArrayType::const_iterator it = points_.begin(); it!=points_.end(); ++it)
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
				//convert input to cgal
				std::vector<Point_2> cgal_points;
				for (typename PointArrayType::const_iterator it = points_.begin(); it!=points_.end(); ++it)
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
				
	      return true;
			}
			
		protected:
			PointArrayType points_;
			typedef CGAL::Cartesian<double>::Point_2 Point_2;
	};
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DCONVEXHULL_H
