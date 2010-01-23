// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

		PeptideIdentifications are grouped by sequence of the respective best-scoring PeptideHit (provided the score is good enough) and retention time data is collected from the "RT" MetaInfo entries.	ID groups with the same sequence in different maps represent points of correspondence between the maps and form the basis of the alignment.
		
		A "consensus retention time scale" is computed from the median retention times (over all maps) of the ID groups. Each map is aligned to this scale as follows:
		The median retention time of each ID group in the map is mapped to the consensus retention time of this group. Cubic spline smoothing is used to convert this mapping to a smooth function. Retention times in the map are transformed to the consensus scale by applying this function.

		@experimental This algorithm is work in progress and has not been thoroughly tested yet.

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
		virtual void alignPeptideIdentifications(
			std::vector<std::vector<PeptideIdentification> >&,
			std::vector<TransformationDescription>&);
		
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

		/// Score threshold for peptide hits
		DoubleReal score_threshold_;
		
		/**
			 @brief Compute the median of a list of values
			 
			 @param values input values (will be sorted)
			 @param sorted are values already sorted? (sorting step can be saved)
			 
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

			 @param peptides input peptide IDs (lists of peptide hits will be sorted)
			 @param rt_data lists of RT values for diff. peptide sequences (output)
		*/
		void getRetentionTimes_(std::vector<PeptideIdentification>& peptides,
														SeqToList& rt_data);

		/**
			 @brief Collect retention time data ("RT" MetaInfo) from peptide IDs annotated to spectra

			 @param experiment Input map for RT data
			 @param rt_data lists of RT values for diff. peptide sequences (output)
		*/
		void getRetentionTimes_(MSExperiment<>& experiment, SeqToList& rt_data);

		/**
			 @brief Collect retention time data ("RT" MetaInfo) from peptide IDs annotated to features

			 Depending on global parameter @p use_unassigned_peptides, unassigned peptide IDs may be used in addition to assigned IDs

			 @param features Input features for RT data
			 @param rt_data lists of RT values for diff. peptide sequences (output)
		*/
		void getRetentionTimes_(FeatureMap<>& features, SeqToList& rt_data);

		/**
			 @brief Compute retention time transformations from RT data grouped by peptide sequence

			 @param rt_data Lists of RT values for diff. peptide sequences, per dataset (input, will be sorted)
			 @param transforms Resulting transformations, per dataset (output)
			 @param sorted Are RT lists already sorted? (see @p median_)
		*/
		void computeTransformations_(std::vector<SeqToList>& rt_data,
																 std::vector<TransformationDescription>&
																 transforms, bool sorted = false);
			
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
