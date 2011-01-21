// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_CONSENSUSMAP_H
#define OPENMS_KERNEL_CONSENSUSMAP_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>

namespace OpenMS
{
	/**
		@brief A container for consensus elements.

		A %ConsensusMap is a container holding 2-dimensional consensus elements (ConsensusFeature)
		which in turn represent combined elements of 2-dimensional experiments.
		The map is implemented as a vector of elements.

		The map indices used in the consensus features should be registered in this class.

		@ingroup Kernel
  */
	class OPENMS_DLLAPI ConsensusMap
		: public std::vector<ConsensusFeature>,
			public MetaInfoInterface,
			public RangeManager<2>,
			public DocumentIdentifier,
			public UniqueIdInterface,
			public UniqueIdIndexer<ConsensusMap>
	{

		public:

			/// Source file description for input files
			struct FileDescription
				: public MetaInfoInterface
			{
				///Default constructor
				FileDescription()
					: MetaInfoInterface(),
						filename(),
						label(),
						size(0),
						unique_id(UniqueIdInterface::INVALID)
				{
				}

				/// File name of the file
				String filename;
				/// Label e.g. 'heavy' and 'light' for ICAT, or 'sample1' and 'sample2' for label-free quantitation
				String label;
				/// @brief Number of elements (features, peaks, ...).
				/// This is e.g. used to check for correct element indices when writing a consensus map TODO fix that
				Size size;
				/// Unique id of the file
				UInt64 unique_id;
			};

			///@name Type definitions
			//@{
			typedef std::vector<ConsensusFeature > Base;
			typedef RangeManager<2> RangeManagerType;
			typedef Map<UInt64,FileDescription> FileDescriptions;
			/// Mutable iterator
			typedef std::vector<ConsensusFeature>::iterator Iterator;
			/// Non-mutable iterator
			typedef std::vector<ConsensusFeature>::const_iterator ConstIterator;
			/// Mutable reverse iterator
			typedef std::vector<ConsensusFeature>::reverse_iterator ReverseIterator;
			/// Non-mutable reverse iterator
			typedef std::vector<ConsensusFeature>::const_reverse_iterator ConstReverseIterator;
			//@}

			/// Default onstructor
			inline ConsensusMap()
				: Base(),
					MetaInfoInterface(),
					RangeManagerType(),
					DocumentIdentifier(),
					UniqueIdInterface(),
					UniqueIdIndexer<ConsensusMap>(),
					file_description_(),
					experiment_type_(),
					protein_identifications_(),
					unassigned_peptide_identifications_(),
					data_processing_()
			{
			}

			/// Copy constructor
			inline ConsensusMap(const ConsensusMap& source)
				: Base(source),
					MetaInfoInterface(source),
					RangeManagerType(source),
					DocumentIdentifier(source),
					UniqueIdInterface(source),
					UniqueIdIndexer<ConsensusMap>(source),
					file_description_(source.file_description_),
					experiment_type_(source.experiment_type_),
					protein_identifications_(source.protein_identifications_),
					unassigned_peptide_identifications_(source.unassigned_peptide_identifications_),
					data_processing_(source.data_processing_)
			{
			}

			/// Destructor
			inline ~ConsensusMap()
			{
			}

			/// Creates a ConsensusMap with n elements
			inline ConsensusMap(Base::size_type n)
				: Base(n),
					MetaInfoInterface(),
					RangeManagerType(),
					DocumentIdentifier(),
					UniqueIdInterface(),
					file_description_(),
					experiment_type_(),
					protein_identifications_(),
					unassigned_peptide_identifications_(),
					data_processing_()
			{
			}

			/// Assignment operator
			ConsensusMap& operator=(const ConsensusMap& source)
			{
				if (this==&source) return *this;

				Base::operator=(source);
				MetaInfoInterface::operator=(source);
				RangeManagerType::operator=(source);
				DocumentIdentifier::operator=(source);
        UniqueIdInterface::operator=(source);
				file_description_ = source.file_description_;
				experiment_type_ = source.experiment_type_;
				protein_identifications_ = source.protein_identifications_;
				unassigned_peptide_identifications_ = source.unassigned_peptide_identifications_;
				data_processing_ = source.data_processing_;

				return *this;
			}

			/**
				@brief Clears all data and meta data
				
				@param clear_meta_data If @em true, all meta data is cleared in addition to the data.
			*/ 
			void clear(bool clear_meta_data = true)
			{
				Base::clear();
					
				if (clear_meta_data)
				{
					clearMetaInfo();
					clearRanges();
					this->DocumentIdentifier::operator=(DocumentIdentifier()); // no "clear" method
					clearUniqueId();
					file_description_.clear();
					experiment_type_.clear();
					protein_identifications_.clear();
					unassigned_peptide_identifications_.clear();
					data_processing_.clear();
				}
			}		

			/// Non-mutable access to the file descriptions
			inline const FileDescriptions& getFileDescriptions() const
			{
				return file_description_;
			}

			/// Mutable access to the file descriptions
			inline FileDescriptions& getFileDescriptions()
			{
				return file_description_;
			}

			/// Non-mutable access to the experiment type
			inline const String& getExperimentType() const
			{
				return experiment_type_;
			}

			/// Mutable access to the experiment type
			inline void setExperimentType(const String& experiment_type)
			{
				experiment_type_ = experiment_type;
			}

#if 0

      // OBSOLETE and BROKEN since we are using unique ids now - the results will be wrong!  (namely, almost always negative)

			/**
				@brief Checks if all map identifiers in FeatureHandles have a filename associated

				@param error_message If the map is not valid, this variable contains the error message
				@return if the map is valid
			*/
			bool isValid(String& error_message) const;
#endif

			/**
				@name Sorting.

				These specialized sorting methods are supported in addition to the standard sorting methods of std::vector. All use stable sorting.
			*/
			//@{
			/// Sorts the peaks according to ascending intensity.
			void sortByIntensity(bool reverse=false)
			{
				if (reverse)
				{
					std::stable_sort(Base::begin(), Base::end(), reverseComparator(ConsensusFeature::IntensityLess()));
				}
				else
				{
					std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::IntensityLess());
				}
			}

			/// Sorts the peaks to RT position.
			void sortByRT()
			{
				std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::RTLess());
			}

			/// Sorts the peaks to m/z position.
			void sortByMZ()
			{
				std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::MZLess());
			}

			/// Lexicographically sorts the peaks by their position (First RT then m/z).
			void sortByPosition()
			{
				std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::PositionLess());
			}

			/// Sorts the peaks according to ascending quality.
			void sortByQuality(bool reverse=false)
			{
				if (reverse)
				{
					std::stable_sort(Base::begin(), Base::end(), reverseComparator(ConsensusFeature::QualityLess()));
				}
				else
				{
					std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::QualityLess());
				}
			}

			/// Sorts with respect to the size (number of elements)
			void sortBySize()
			{
				std::stable_sort(Base::begin(), Base::end(), reverseComparator(ConsensusFeature::SizeLess()));
			}

			/// Sorts with respect to the sets of maps covered by the consensus features (lexicographically).
			void sortByMaps()
			{
				std::stable_sort(Base::begin(), Base::end(), ConsensusFeature::MapsLess());
			}

			//@}

			/**
				@brief Convert a FeatureMap (of any feature type) to a ConsensusMap.  

				Each ConsensusFeature contains a map index, so this has to be given as
				well. The previous content of @p output_map is cleared. An arguable
				design decision is that the unique id of the FeatureMap is copied (!) to
				the ConsensusMap, because that is the way it is meant to be used in the
				algorithms.

				Only the first (!) @p n elements are copied. (This parameter exists
				mainly for compatibility with @p convert for MSExperiments. To use it in
				a meaningful way, apply one of the sorting methods to @p input_map
				beforehand.)
 
				@param input_map_index The index of the input map.
				@param input_map The container to be converted.
				@param output_map The resulting ConsensusMap.
				@param n The maximum number of elements to be copied.
			*/
			template <typename FeatureT>
				static void convert(UInt64 const input_map_index, 
														FeatureMap<FeatureT> const & input_map, 
														ConsensusMap& output_map, Size n = -1)
      {
				if ( n > input_map.size() )
				{
					n = input_map.size();
				}

				output_map.clear(true);
				output_map.reserve(n);
				
				// An arguable design decision, see above.
        output_map.setUniqueId(input_map.getUniqueId());

        for (UInt64 element_index = 0; element_index < n; ++element_index )
				{
					output_map.push_back( ConsensusFeature( input_map_index, input_map[element_index] ) );
				}
				output_map.getFileDescriptions()[input_map_index].size = (Size) input_map.size();
        output_map.setProteinIdentifications(input_map.getProteinIdentifications());
        output_map.setUnassignedPeptideIdentifications(input_map.getUnassignedPeptideIdentifications());
				output_map.updateRanges();
			}

			/**
				@brief Similar to @p convert for FeatureMaps.

				Only the @p n most intense elements are copied.

        Currently MSExperiment<> does not have a unique id but ConsensusMap has
        one, so we assign a new one here.

				@param input_map_index The index of the input map.
				@param input_map The input map to be converted.
				@param output_map The resulting ConsensusMap.
				@param n The maximum number of elements to be copied.
			*/
			static void convert(UInt64 const input_map_index, MSExperiment<> & input_map, ConsensusMap& output_map, Size n = -1)
			{
				output_map.clear(true);

				// see @todo above
        output_map.setUniqueId();

				input_map.updateRanges(1);
				if ( n > input_map.getSize() )
				{
					n = input_map.getSize();
				}
				output_map.reserve(n);
				std::vector<Peak2D> tmp;
				tmp.reserve(input_map.getSize());
				input_map.get2DData(tmp); // TODO Avoid tripling the memory consumption by this call
				std::partial_sort( tmp.begin(), tmp.begin()+n, tmp.end(), reverseComparator(Peak2D::IntensityLess()) );
				for (Size element_index = 0; element_index < n; ++element_index )
				{
					output_map.push_back( ConsensusFeature(input_map_index, tmp[element_index], element_index ) );
				}
				output_map.getFileDescriptions()[input_map_index].size = n;
				output_map.updateRanges();
			}

      /**
				@brief Convert a ConsensusMap to a FeatureMap (of any feature type).

				The previous content of output_map is cleared. UID's of the elements and
				the container is copied if the @p keep_uids flag is set.

				@param input_map The container to be converted.
        @param keep_uids Shall the UID's of the elements and the container be kept or created anew
				@param output_map The resulting ConsensusMap.
			*/
      template <typename FeatureT>
      static void convert(ConsensusMap const& input_map, const bool keep_uids, FeatureMap<FeatureT> & output_map )
      {
        output_map.clear(true);
				output_map.resize(input_map.size());
				// fm.MetaInfoInterface::operator=(input_map); // not available ...
				output_map.DocumentIdentifier::operator=(input_map);
				if (keep_uids) output_map.UniqueIdInterface::operator=(input_map);
        else output_map.setUniqueId();
				output_map.setProteinIdentifications(input_map.getProteinIdentifications());
				output_map.setUnassignedPeptideIdentifications(input_map.getUnassignedPeptideIdentifications());
				for ( Size i = 0; i < input_map.size(); ++i )
				{
					Feature & f = output_map[i];
					const ConsensusFeature & c = input_map[i];
					f.BaseFeature::operator=(c);
          if (!keep_uids) f.setUniqueId();
				}

      }

			// Docu in base class
			void updateRanges();

			/// Swaps the content of this map with the content of @p from
			void swap(ConsensusMap& from)
			{
				ConsensusMap tmp;

				//swap range information
				tmp.RangeManagerType::operator=(*this);
				this->RangeManagerType::operator=(from);
				from.RangeManagerType::operator=(tmp);

				//swap consensus features
				Base::swap(from);

				// swap DocumentIdentifier
				DocumentIdentifier::swap(from);

				// swap unique id
				UniqueIdInterface::swap(from);

        // swap unique id index
        UniqueIdIndexer<ConsensusMap>::swap(from);

				// swap the remaining members
				std::swap(file_description_, from.file_description_);
				experiment_type_.swap(from.experiment_type_);
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

			/// Equality operator
			bool operator == (const ConsensusMap& rhs) const
			{
				return
					std::operator==(*this, rhs) &&
					MetaInfoInterface::operator==(rhs) &&
					RangeManagerType::operator==(rhs) &&
					DocumentIdentifier::operator==(rhs) &&
					UniqueIdInterface::operator == (rhs) &&
					file_description_ == rhs.file_description_ &&
					experiment_type_ == rhs.experiment_type_ &&
					protein_identifications_==rhs.protein_identifications_ &&
					unassigned_peptide_identifications_==rhs.unassigned_peptide_identifications_ &&
					data_processing_ == rhs.data_processing_
					;
			}

			/// Equality operator
			bool operator != (const ConsensusMap& rhs) const
			{
				return !(operator==(rhs));
			}

      /**@brief Applies a member function of Type to the container itself and all consensus features.
         The returned values are accumulated.

         <b>Example:</b>  The following will print the number of features with invalid unique ids (plus 1 if the container has an invalid UID as well):
         @code
         ConsensusMap cm;
         (...)
         std::cout << cm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId) << std::endl;
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
          assignments += ((*iter).*member_function)();
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
          assignments += ((*iter).*member_function)();
        }
        return assignments;
      }

		protected:

			/// Map from index to file description
			FileDescriptions file_description_;

			/// type of experiment (label-free, itraq, ...); see xsd schema
			String experiment_type_;

			/// protein identifications
			std::vector<ProteinIdentification> protein_identifications_;

			/// protein identifications
			std::vector<PeptideIdentification> unassigned_peptide_identifications_;

			/// applied data processing
			std::vector<DataProcessing> data_processing_;

	};

  ///Print the contents of a ConsensusMap to a stream.
  OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const ConsensusMap& cons_map);

} // namespace OpenMS

#endif // OPENMS_KERNEL_CONSENSUSMAP_H
