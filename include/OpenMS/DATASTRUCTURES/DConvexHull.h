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
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

#include <vector>

namespace OpenMS
{
	/**	
		@brief A d-dimensional convex hull representation
		
		@todo improve operator== and addPoint (Marc)
		
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
				for (typename PointArrayType::const_iterator it = points.begin(); it!=points.end(); ++it)
				{
					addPoint(*it);
				}
				
				return *this;
			}

			/// equality operator
			bool operator==(const DConvexHull& rhs) const
			{
				//TODO check if all point are contained in the other chonvex hull
				return points_ == rhs.points_;
			}
			
			/// adds a point
			void addPoint(const PointType& point)
			{
				//TODO check if already contained
				points_.push_back(point);
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
			
			/// returns if the @p point lies in the convex hull
			bool encloses(const PointType& point) const
			{
				//TODO
				return false;
			}
			
		protected:
			PointArrayType points_;
	};
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DCONVEXHULL_H
