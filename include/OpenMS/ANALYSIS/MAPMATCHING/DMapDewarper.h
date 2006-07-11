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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_DMAPDEWARPER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DMAPDEWARPER_H

#include<OpenMS/ANALYSIS/MAPMATCHING/DBaseMapping.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>

#include<OpenMS/KERNEL/DFeatureMap.h>

#include<vector>

namespace OpenMS
{
	
	/**
	 	@brief This class applies a transformation as computed by
	 	class DBaseMapMatcher to a set of features.
	 		 	
	 **/
	template <UnsignedInt D, typename Traits = KernelTraits>		
	class DMapDewarper
	{
		public:
			
		/** @name Type definitions
		*/
		//@{	
		/// The grid is simply a vector of cells.
		typedef DGrid<D> Grid;
		/// 
		typedef DFeatureMap<D> FeatureMapType;
		///
		typedef DBaseMapping<1> MappingType;
		///
		typedef std::vector<MappingType*> MappingVector;
		//@}
		
		/// Constructor
		DMapDewarper()
		 : grid_(), features_() {}
		
		/// Copy constructor
		DMapDewarper(const DMapDewarper& source)
			: grid_(source.grid_),
		  	features_(source.features_)  {}
		
		///  Assignment operator
    DMapDewarper& operator = (const DMapDewarper& source)
    {
    	if (&source == this) return *this;
			
			grid_     = source.grid_;
			features_ = source.features_;
			return *this;
		} 
		
		/// equality operator
		bool operator == (const DMapDewarper& rhs)
		{
			return (grid_     == rhs.grid_ &&
							features_ == rhs.features_);
		}	
		
		/// Destructor
		virtual ~DMapDewarper() {}
		
		/* 
		  @brief Dewarps the feature map. 
		  
		  This is a bit ugly. The mapping is D-dimensional function of
		  DPosition<D>. But we work with two one dimensional functions
		  (one for rt and one for m/z). Therefore we take the DPosition<2>
		  of each feature and split it into two instances of DPosition<1> ....
		  
		*/
		void dewarp()
		{
			
			// iterate over all features...
			typename FeatureMapType::iterator feat_iter = features_.begin();
			while (feat_iter != features_.end() )
			{
				// Test in which cell this feature is included
				// and apply the corresponding transformation
				typename Grid::iterator grid_iter = grid_.begin();
				while (grid_iter != grid_.end() )
				{
					if (grid_iter->encloses(feat_iter->getPosition() ) )
					{
						DPosition<D> pos         = feat_iter->getPosition();
						// apply transform for every coordinate
						for (unsigned int i=0; i<D; i++)
						{
							DPosition<1> temp;
							temp[0] = pos[i];
							grid_iter->getMappings().at(i)->apply(temp);
							pos[i] = temp[0];
						}
						feat_iter->setPosition(pos);
					}
					grid_iter++;
				
				} // end while (grid)	
					
				feat_iter++;
			} // end while (features)
			
		}
		
		/** @name Accesssor methods
		*/
		//@{	
		/// Set grid
    void setGrid(Grid& g) { grid_ = g; }
    /// Get grid 
    Grid& getGrid() { return grid_; }
    /// Get grid (non-mutable)
    Grid& getGrid() const { return grid_; }
    /// Set features
    void setFeatures(FeatureMapType& feat) { features_ = feat; }
    /// Get grid 
    FeatureMapType& getFeatures() { return features_; }
    /// Get grid (non-mutable)
    FeatureMapType& getFeatures() const { return features_; }
    //@} 
   
		protected:		
		/// Vector of DRange instances defining a grid over the map
		Grid grid_;
		
		/// The feature that we want to dewarp
		FeatureMapType features_;
					
	}; // end of class DMapDewarper
	
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DMAPDEWARPER_H
