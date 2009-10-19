// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//cl
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
// $Maintainer: Chris Bielow, Clemens Groepl $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_FEATUREMAP_H
#define OPENMS_KERNEL_FEATUREMAP_H

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>

#include <algorithm>
#include <vector>
#include <exception>

namespace OpenMS
{

	/**	
		@brief A container for features.
		
		A map is a container holding 2-dimensional features,
		which in turn represent chemical entities (peptides, proteins, etc.) found
		in a 2-dimensional experiment.
		
		Maps are implemented as vectors of features and have basically the same interface
		as an STL vector has (model of Random Access Container and Back Insertion Sequence).
		
		Feature maps are typically created from peak data of 2D runs through the FeatureFinder.
		
		@ingroup Kernel
	*/
	template <typename FeatureT = Feature >
	class FeatureMap
		: public std::vector<FeatureT>,
			public RangeManager<2>,
			public DocumentIdentifier,
			public UniqueIdInterface,
			public UniqueIdIndexer<FeatureMap<FeatureT> >
	{
		public:
			/**	
				 @name Type definitions
			*/
			//@{
			typedef FeatureT FeatureType;
			typedef RangeManager<2> RangeManagerType;
			typedef std::vector<FeatureType> Base;
			typedef typename Base::iterator Iterator;
			typedef typename Base::const_iterator ConstIterator;
			typedef typename Base::reverse_iterator ReverseIterator;
			typedef typename Base::const_reverse_iterator ConstReverseIterator;
			typedef FeatureType& Reference;
			typedef const FeatureType& ConstReference;
			//@}
			/**	
				 @name Constructors and Destructor
			*/
			//@{
			
			/// Default constructor
			FeatureMap()
				: Base(),
					RangeManagerType(),
					DocumentIdentifier(),
					UniqueIdInterface(),
					UniqueIdIndexer< FeatureMap<FeatureT> >(),
					protein_identifications_(),
					unassigned_peptide_identifications_(),
					data_processing_() 
			{
			}
			
			/// Copy constructor
			FeatureMap(const FeatureMap& source) 
				: Base(source),
					RangeManagerType(source),
					DocumentIdentifier(source),
					UniqueIdInterface(source),
					UniqueIdIndexer< FeatureMap<FeatureT> >(source),
					protein_identifications_(source.protein_identifications_),
					unassigned_peptide_identifications_(source.unassigned_peptide_identifications_),
					data_processing_(source.data_processing_)
			{
			}
			
			/// Destructor
			virtual ~FeatureMap()
			{
			}
			
			//@}
				
			/// Assignment operator
			FeatureMap& operator = (const FeatureMap& rhs)
			{
				if (&rhs==this) return *this;
					
				Base::operator=(rhs);
				RangeManagerType::operator=(rhs);
				DocumentIdentifier::operator=(rhs);
	      UniqueIdInterface::operator = (rhs);
				protein_identifications_ = rhs.protein_identifications_;
				unassigned_peptide_identifications_ = rhs.unassigned_peptide_identifications_;
				data_processing_ = rhs.data_processing_;

				return *this;
			}
	
			/// Equality operator
			bool operator == (const FeatureMap& rhs) const
			{
				return
					std::operator==(*this, rhs) &&
					RangeManagerType::operator==(rhs) &&
					DocumentIdentifier::operator==(rhs) &&
					UniqueIdInterface::operator == (rhs) &&
					protein_identifications_==rhs.protein_identifications_ &&
					unassigned_peptide_identifications_==rhs.unassigned_peptide_identifications_ &&
					data_processing_ == rhs.data_processing_
					;
			}
				
			/// Equality operator
			bool operator != (const FeatureMap& rhs) const
			{
				return !(operator==(rhs));
			}

			/**	
				@name Sorting.
				These simplified sorting methods are supported in addition to	
				the standard sorting methods of std::vector.
			*/
			//@{
			/// Sorts the peaks according to ascending intensity.
			void sortByIntensity(bool reverse=false)
			{ 
				if (reverse)
				{
					std::sort(this->begin(), this->end(), reverseComparator(typename FeatureType::IntensityLess()) );
				}
				else
				{
					std::sort(this->begin(), this->end(), typename FeatureType::IntensityLess() ); 
				}
			}
				
			///Sort features by position. Lexicographical comparison (first RT then m/z) is done.
			void sortByPosition() 
			{ 
				std::sort(this->begin(), this->end(), typename FeatureType::PositionLess() );
			}
			
			///Sort features by RT position.
			void sortByRT() 
			{ 
				std::sort(this->begin(), this->end(), typename FeatureType::RTLess() );
			}

			///Sort features by m/z position.
			void sortByMZ() 
			{ 
				std::sort(this->begin(), this->end(), typename FeatureType::MZLess() );
			}
			
			///Sort features by ascending overall quality.
			void sortByOverallQuality(bool reverse=false) 
			{
				if (reverse)
				{
					std::sort(this->begin(), this->end(), reverseComparator(typename FeatureType::OverallQualityLess()) );
				}
				else
				{
					std::sort(this->begin(), this->end(), typename FeatureType::OverallQualityLess() );
				}
			}
			//@}
			
			// Docu in base class
			void updateRanges()
			{
				this->clearRanges();
				updateRanges_(this->begin(),this->end());
				
				//enlarge the range by the convex hull points
				for (Size i=0; i<this->size(); ++i)
				{
					DBoundingBox<2> box = this->operator[](i).getConvexHull().getBoundingBox();
					if (!box.isEmpty())
					{
						//update RT
						if (box.min()[Peak2D::RT] < this->pos_range_.min()[Peak2D::RT])
						{
							this->pos_range_.setMinX(box.min()[Peak2D::RT]);
						}
						if (box.max()[Peak2D::RT] > this->pos_range_.max()[Peak2D::RT])
						{
							this->pos_range_.setMaxX(box.max()[Peak2D::RT]);
						}
						//update m/z
						if (box.min()[Peak2D::MZ] < this->pos_range_.min()[Peak2D::MZ])
						{
							this->pos_range_.setMinY(box.min()[Peak2D::MZ]);
						}
						if (box.max()[Peak2D::MZ] > this->pos_range_.max()[Peak2D::MZ])
						{
							this->pos_range_.setMaxY(box.max()[Peak2D::MZ]);
						}
					}
				}
			}

			/// Swaps the content of this map with the content of @p from
			void swap(FeatureMap& from)
			{
				FeatureMap tmp;

        // swap the actual features
        Base::swap(from);

				// swap range information
				tmp.RangeManagerType::operator=(*this);
				this->RangeManagerType::operator=(from);
				from.RangeManagerType::operator=(tmp);

				// swap DocumentIdentifier
				DocumentIdentifier::swap(from);

        // swap unique id
        UniqueIdInterface::swap(from);

        // swap unique id index
        UniqueIdIndexer<FeatureMap<FeatureT> >::swap(from);
				
				// swap the remaining members
				protein_identifications_.swap(from.protein_identifications_);
				unassigned_peptide_identifications_.swap(from.unassigned_peptide_identifications_);
				data_processing_.swap(from.data_processing_);
			}
			
			/// non-mutable access to the protein identifications
		 	const std::vector<ProteinIdentification>& getProteinIdentifications() const
		 	{
		  	return protein_identifications_;	   		
		 	}	
		 		    	
			/// mutable access to the protein identifications
		  std::vector<ProteinIdentification>& getProteinIdentifications()
		  {
		  	return protein_identifications_;	
		  }

			/// sets the protein identifications
		  void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications)
		  {
		  	protein_identifications_ = protein_identifications;
		  }
		  
			/// non-mutable access to the unassigned peptide identifications
			const std::vector<PeptideIdentification>& getUnassignedPeptideIdentifications() const
			{
				return unassigned_peptide_identifications_;	   		
			}	
			
			/// mutable access to the unassigned peptide identifications
			std::vector<PeptideIdentification>& getUnassignedPeptideIdentifications()
			{
				return unassigned_peptide_identifications_;	
			}
			
			/// sets the unassigned peptide identifications
			void setUnassignedPeptideIdentifications(const std::vector<PeptideIdentification>& unassigned_peptide_identifications)
			{
				unassigned_peptide_identifications_ = unassigned_peptide_identifications;
			}

			/// returns a const reference to the description of the applied data processing 
			const std::vector<DataProcessing>& getDataProcessing() const
			{
				return data_processing_; 
			}

			/// returns a mutable reference to the description of the applied data processing 
			std::vector<DataProcessing>& getDataProcessing()
			{
				return data_processing_; 
			}
			
			/// sets the description of the applied data processing 
			void setDataProcessing(const std::vector<DataProcessing>& processing_method)
			{
				data_processing_ = processing_method; 
			}

			/**
				@brief Clears all data and meta data
				
				@param clear_meta_data If @em true, all meta data is cleared in addition to the data.
			*/ 
			void clear(bool clear_meta_data)
			{
				Base::clear();
				
				if (clear_meta_data)
				{
					clearRanges();
					this->DocumentIdentifier::operator=(DocumentIdentifier()); // no "clear" method
					clearUniqueId();
					protein_identifications_.clear();
					unassigned_peptide_identifications_.clear();
					data_processing_.clear();
				}
			}
		
      /**@brief Applies a member function of Type to all features, including subordinates.
         The returned values are accumulated.

         <b>Example:</b>  The following will print the number of features with invalid unique ids:
         @code
         FeatureMap<> fm;
         (...)
         std::cout << fm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId) << std::endl;
         @endcode
         See e.g. UniqueIdInterface for what else can be done this way.
      */
      template < typename Type >
      Size applyMemberFunction( Size (Type::*member_function)() )
      {
        Size assignments = 0;
        assignments += ((*this).*member_function)();
        for ( Iterator iter = this->begin(); iter != this->end(); ++iter)
        {
          assignments += iter->applyMemberFunction(member_function);
        }
        return assignments;
      }

      /// The "const" variant.
      template < typename Type >
      Size applyMemberFunction( Size (Type::*member_function)() const ) const
      {
        Size assignments = 0;
        assignments += ((*this).*member_function)();
        for ( ConstIterator iter = this->begin(); iter != this->end(); ++iter)
        {
          assignments += iter->applyMemberFunction(member_function);
        }
        return assignments;
      }

		protected:

			/// protein identifications
			std::vector<ProteinIdentification> protein_identifications_;

			/// protein identifications
			std::vector<PeptideIdentification> unassigned_peptide_identifications_;
			
			/// applied data processing
			std::vector<DataProcessing> data_processing_;
	};
	
	/// Print content of a feature map to a stream.
	template <typename FeatureType >
	std::ostream& operator << (std::ostream& os, const FeatureMap<FeatureType>& map)
	{
		os << "# -- DFEATUREMAP BEGIN --"<< std::endl;
		os << "# POS \tINTENS\tOVALLQ\tCHARGE\tUniqueID" << std::endl;
		for (typename FeatureMap<FeatureType>::const_iterator iter = map.begin(); iter!=map.end(); iter++)
		{
			os << iter->getPosition() << '\t'
				 << iter->getIntensity() << '\t'
				 << iter->getOverallQuality() << '\t'
				 << iter->getCharge() << '\t'
				 << iter->getUniqueId()
				 << std::endl;
		}
		os << "# -- DFEATUREMAP END --"<< std::endl;
		return os;
	}
	
} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATUREMAP_H
