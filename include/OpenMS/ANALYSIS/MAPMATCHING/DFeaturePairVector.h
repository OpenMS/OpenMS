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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_DFEATUREPAIRVECTOR_H
#define OPENMS_ANALYSIS_MAPMATCHING_DFEATUREPAIRVECTOR_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>

#include <vector>

namespace OpenMS
{

	/** @brief Data structure containing pairs of features (as e.g. produced by DMapMatcherRegression)
	 */	
	template <Size D, typename FeatureT = DFeature<D> >
	class DFeaturePairVector
		: public std::vector< DFeaturePair<D, FeatureT > >
	{
	 public:
			
		/**	
			 @name Type definitions
		*/
		//@{
		typedef DFeaturePair<D,FeatureT> PairType;
		typedef std::vector<PairType> Base;
		typedef std::vector<PairType> ContainerType;
		typedef typename ContainerType::iterator Iterator;
		typedef typename ContainerType::const_iterator ConstIterator;
		typedef typename ContainerType::reverse_iterator ReverseIterator;
		typedef typename ContainerType::const_reverse_iterator ConstReverseIterator;
		typedef PairType & Reference;
		typedef const PairType & ConstReference;
			
		// STL compatibility
		typedef PairType value_type;
		typedef PairType* pointer;
		typedef const PairType* const_pointer;
		typedef Reference reference;
		typedef ConstReference const_reference;
		typedef typename ContainerType::size_type size_type;
		typedef typename ContainerType::difference_type difference_type;
		typedef Iterator iterator;
		typedef ConstIterator const_iterator;
		typedef ReverseIterator reverse_iterator;
		typedef ConstReverseIterator const_reverse_iterator;
		//@}
			
		/**	
			 @name Constructors and Destructor
		*/
		//@{
		/// Default constructor
		DFeaturePairVector() {}
		/// Copy constructor
		DFeaturePairVector(const DFeaturePairVector& vec) :
			Base(vec)
		{}
		/// Destructor
		virtual ~DFeaturePairVector() {}
		//@}
			
		/// Assignment operator
		DFeaturePairVector& operator = (const DFeaturePairVector& rhs)
		{
			if (&rhs==this) return *this;
				
			Base::operator=(rhs);
								
			return *this;
		}

		/// Equality operator
		bool operator == (const DFeaturePairVector& rhs) const
		{
			return	std::operator==(*this, rhs)	;				
		}
			
		/// Unequality operator
		bool operator != (const DFeaturePairVector& rhs) const
		{
			return !(operator==(rhs));
		}
						
	}; // end of class DFeaturePairVector
	
	///Print the contents to a stream.
	template <Size D>
	std::ostream& operator << (std::ostream& os, const DFeaturePairVector<D>& pairs)
	{
		os << "-- DFEATUREPAIRVECTOR BEGIN --" << std::endl;
		for (typename DFeaturePairVector<D>::const_iterator it = pairs.begin(); it!=pairs.end(); ++it)
		{
			os << "FIRST: " << std::endl;
			os << it->getFirst() << std::endl;
			
			os << "SECOND: " << std::endl;
			os << it->getSecond() << std::endl;
		}
		os << "-- DFEATUREPAIRVECTOR END --"<<std::endl;
		return os;
	}
	
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DFEATUREPAIRVECTOR_H
