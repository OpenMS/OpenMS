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
// $Maintainer: Hendrik Weisser $
// $Authors: Eva Lange, Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMIDENTIFICATION_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMIDENTIFICATION_H

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <map>

namespace OpenMS
{
	/**
		@brief A map alignment algorithm based on peptide identifications from MS2 spectra.

		PeptideIdentification instances are grouped by sequence of the respective best-scoring PeptideHit (provided the score is good enough) and retention time data is collected from the "RT" MetaInfo entries.	ID groups with the same sequence in different maps represent points of correspondence between the maps and form the basis of the alignment.
		
		Each map is aligned to a reference retention time scale. This time scale can either come from a reference file (@p reference parameter) or be computed as a consensus of the input maps (median retention times over all maps of the ID groups). The maps are then aligned to this scale as follows:\n
		The median retention time of each ID group in a map is mapped to the reference retention time of this group. Cubic spline smoothing is used to convert this mapping to a smooth function. Retention times in the map are transformed to the consensus scale by applying this function.

	  @htmlinclude OpenMS_MapAlignmentAlgorithmIdentification.parameters

		@ingroup MapAlignment
	*/
	class OPENMS_DLLAPI MapAlignmentAlgorithmIdentification
		: public MapAlignmentAlgorithm
	{
	 public:
		/// Default constructor
		MapAlignmentAlgorithmIdentification();

		/// Destructor
		virtual ~MapAlignmentAlgorithmIdentification();

		// Docu in base class
		virtual void alignPeakMaps(std::vector<MSExperiment<> >&,
															 std::vector<TransformationDescription>&);

		// Docu in base class
		virtual void alignFeatureMaps(std::vector<FeatureMap<> >&,
																	std::vector<TransformationDescription>&);

		// Docu in base class
		virtual void alignConsensusMaps(std::vector<ConsensusMap>&,
																		std::vector<TransformationDescription>&);

		// Docu in base class
		virtual void alignPeptideIdentifications(
			std::vector<std::vector<PeptideIdentification> >&,
			std::vector<TransformationDescription>&);

		// Docu in base class
		virtual void setReference(Size reference_index=0, 
															const String& reference_file="");

		// Docu in base class
		virtual void getDefaultModel(String& model_type, Param& params);
	
		/// Creates a new instance of this class (for Factory)
		static MapAlignmentAlgorithm* create()
		{
			return new MapAlignmentAlgorithmIdentification();
		}

		/// Returns the product name (for the Factory)
		static String getProductName()
		{
			return "identification";
		}

	 protected:

		/// Type to store retention times given for individual peptide sequences
		typedef std::map<String, DoubleList> SeqToList;
		
		/// Type to store one representative retention time per peptide sequence
		typedef std::map<String, DoubleReal> SeqToValue;

		/// Index of input file to use as reference (1-based!)
		Size reference_index_;

		/// Reference retention times (per peptide sequence)
		SeqToValue reference_;

		/// Score threshold for peptide hits
		DoubleReal score_threshold_;
		
		/**
			 @brief Compute the median of a list of values
			 
			 @param values Input values (will be sorted)
			 @param sorted Are values already sorted? (sorting step can be saved)
			 
			 @throw Exception::IllegalArgument if the input list is empty
		*/
		DoubleReal median_(DoubleList& values, bool sorted = false);

		/**
			 @brief Compute the median retention time for each peptide sequence
			 
			 @param rt_data Lists of RT values for diff. peptide sequences (input, will be sorted)
			 @param medians Median RT values for the peptide sequences (output)
			 @param sorted Are RT lists already sorted? (see @p median_)
			 
			 @throw Exception::IllegalArgument if the input list is empty
		*/
		void computeMedians_(SeqToList& rt_data, SeqToValue& medians,
												 bool sorted = false);

		/// Check if peptide ID contains a hit that passes the significance threshold @p score_threshold_ (list of peptide hits will be sorted)
		bool hasGoodHit_(PeptideIdentification& peptide);

		/**
			 @brief Collect retention time data ("RT" MetaInfo) from peptide IDs

			 @param peptides Input peptide IDs (lists of peptide hits will be sorted)
			 @param rt_data Lists of RT values for diff. peptide sequences (output)
		*/
		void getRetentionTimes_(std::vector<PeptideIdentification>& peptides,
														SeqToList& rt_data);

		/**
			 @brief Collect retention time data ("RT" MetaInfo) from peptide IDs annotated to spectra

			 @param experiment Input map for RT data
			 @param rt_data Lists of RT values for diff. peptide sequences (output)
		*/
		void getRetentionTimes_(MSExperiment<>& experiment, SeqToList& rt_data);

		/**
			 @brief Collect retention time data ("RT" MetaInfo) from peptide IDs contained in feature maps or consensus maps

			 The following global flags (mutually exclusive) influence the processing:\n
			 Depending on @p use_unassigned_peptides, unassigned peptide IDs are used in addition to IDs annotated to features.\n
			 Depending on @p use_feature_rt, feature retention times are used instead of peptide retention times.

			 @param features Input features for RT data
			 @param unassigned List of unassigned peptide IDs
			 @param rt_data Lists of RT values for diff. peptide sequences (output)
		*/
		template <typename MapType> void getRetentionTimes_(MapType& features,
																												SeqToList& rt_data);

		/**
			 @brief Compute retention time transformations from RT data grouped by peptide sequence

			 @param rt_data Lists of RT values for diff. peptide sequences, per dataset (input, will be sorted)
			 @param transforms Resulting transformations, per dataset (output)
			 @param sorted Are RT lists already sorted? (see @p median_)
		*/
		void computeTransformations_(std::vector<SeqToList>& rt_data,
																 std::vector<TransformationDescription>&
																 transforms, bool sorted = false);

		/**
			 @brief Check that parameter values are valid

			 Currently only 'min_run_occur' is checked.

			 @param runs Number of runs (input files) to be aligned
		*/
		void checkParameters_(const Size runs);
			
		/**
			 @brief Get reference retention times

			 If a reference file is supplied via the @p reference parameter, extract retention time information and store it in #reference_.
		*/
		void getReference_();

		/**
			 @brief Align feature maps or consensus maps

			 Helper function containing common functionality for alignFeatureMaps and alignConsensusMaps.
		*/
		template <typename MapType> void alignFeatureOrConsensusMaps_(
			std::vector<MapType>& maps, std::vector<TransformationDescription>& 
			transformations);

	 private:

		/// Copy constructor intentionally not implemented -> private
		MapAlignmentAlgorithmIdentification(
			const MapAlignmentAlgorithmIdentification&);

		///Assignment operator intentionally not implemented -> private
		MapAlignmentAlgorithmIdentification& operator=(
			const MapAlignmentAlgorithmIdentification&);

	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMIDENTIFICATION_H
