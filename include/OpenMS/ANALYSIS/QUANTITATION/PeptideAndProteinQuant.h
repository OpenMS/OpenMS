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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_PEPTIDEANDPROTEINQUANT_H
#define OPENMS_ANALYSIS_QUANTITATION_PEPTIDEANDPROTEINQUANT_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

namespace OpenMS
{
	/**
		@brief Helper class for peptide and protein quantification based on feature data annotated with IDs

		This class is used by @ref TOPP_ProteinQuantifier. See there for further documentation.

		@htmlinclude OpenMS_PeptideAndProteinQuant.parameters
	*/
	class OPENMS_DLLAPI PeptideAndProteinQuant: public DefaultParamHandler
	{
	public:
		
		/// Mapping: sample ID -> abundance
		typedef std::map<UInt64, DoubleReal> SampleAbundances;

		/// Quantitative and associated data for a peptide
		struct PeptideData
		{
			/// mapping: charge -> sample -> abundance
			std::map<Int, SampleAbundances> abundances;

			/// mapping: sample -> total abundance
			SampleAbundances total_abundances;

			/// protein accessions for this peptide
			std::set<String> accessions;

			/// number of identifications
			Size id_count;

			/// constructor
			PeptideData(): id_count(0) {}
		};
		
		/// Mapping: peptide sequence (modified) -> peptide data
		typedef std::map<AASequence, PeptideData> PeptideQuant;

		/// Quantitative and associated data for a protein
		struct ProteinData
		{
			/// mapping: peptide (unmodified) -> sample -> abundance
			std::map<String, SampleAbundances> abundances;

			/// mapping: sample -> total abundance
			SampleAbundances total_abundances;

			/// total number of identifications (of peptides mapping to this protein)
			Size id_count;

			/// constructor
			ProteinData(): id_count(0) {}
		};

		/// Mapping: protein accession -> protein data
		typedef std::map<String, ProteinData> ProteinQuant;

		/// Statistics for processing summary
		struct Statistics 
		{
			/// number of samples
			Size n_samples;

			/// protein statistics
			Size quant_proteins, too_few_peptides;

			/// peptide statistics
			Size quant_peptides, total_peptides;

			/// feature statistics
			Size quant_features, total_features, blank_features, ambig_features;

			/// constructor
			Statistics(): n_samples(0), quant_proteins(0), too_few_peptides(0), 
										quant_peptides(0), total_peptides(0), quant_features(0), 
										total_features(0), blank_features(0), ambig_features(0) {}
		}; 

		/// Constructor
		PeptideAndProteinQuant();

		/// Destructor
		~PeptideAndProteinQuant() {};

		/**
			 @brief Compute peptide abundances from data in a feature map
			 
			 Parameters should be set before using this method, as setting parameters will clear all results.
		*/
		void quantifyPeptides(FeatureMap<>& features);

		/**
			 @brief Compute peptide abundances from data in a consensus map
			 
			 Parameters should be set before using this method, as setting parameters will clear all results.
		*/
		void quantifyPeptides(ConsensusMap& consensus);

		/**
			 @brief Compute protein abundances.

			 Peptide abundances must be computed first with @p quantifyPeptides. Optional information about groups of indistinguishable proteins (from ProteinProphet) can be supplied via @p proteins.
		*/
		void quantifyProteins(const ProteinIdentification& proteins = 
													ProteinIdentification());

		/// Get summary statistics
		const Statistics& getStatistics();

		/// Get peptide abundance data
		const PeptideQuant& getPeptideResults();

		/// Get protein abundance data
		const ProteinQuant& getProteinResults();

	private:

		/// Processing statistics for output in the end
		Statistics stats_;

		/// Peptide quantification data
		PeptideQuant pep_quant_;

		/// Protein quantification data
		ProteinQuant prot_quant_;


		/**
			 @brief Get the "canonical" annotation (a single peptide hit) of a feature/consensus feature from the associated list of peptide identifications.

			 Only the best-scoring peptide hit of each ID in @p peptides is taken into account. The hits of each ID must already be sorted! If there's more than one ID and the best hits are not identical by sequence, or if there's no peptide ID, an empty peptide hit (for "ambiguous/no annotation") is returned.
			 Protein accessions from identical peptide hits are accumulated.
		*/
		PeptideHit getAnnotation_(std::vector<PeptideIdentification>& peptides);

		/**
			 @brief Gather quantitative information from a feature.

			 Store quantitative information from @p feature in member @p pep_quant_, based on the peptide annotation in @p hit. If @p hit is empty ("ambiguous/no annotation"), nothing is stored.
		*/
		void quantifyFeature_(const FeatureHandle& feature, const PeptideHit& hit);

		/**
			 @brief Order keys (charges/peptides for peptide/protein quantification) according to how many samples they allow to quantify, breaking ties by total abundance.

			 The keys of @p abundances are stored ordered in @p result, best first.
		*/
		template <typename T>
		void orderBest_(const std::map<T, SampleAbundances> abundances, 
										std::vector<T>& result)
		{
			typedef std::pair<Size, DoubleReal> PairType;
			std::multimap<PairType, T, std::greater<PairType> > order;
			for (typename std::map<T, SampleAbundances>::const_iterator ab_it = 
						 abundances.begin(); ab_it != abundances.end(); ++ab_it)
			{
				DoubleReal total = 0.0;
				for (SampleAbundances::const_iterator samp_it = ab_it->second.begin();
						 samp_it != ab_it->second.end(); ++samp_it)
				{
					total += samp_it->second;
				}
				if (total <= 0.0) continue; // not quantified
				PairType key = std::make_pair(ab_it->second.size(), total);
				order.insert(std::make_pair(key, ab_it->first));
			}
			result.clear();
			for (typename std::multimap<PairType, T, std::greater<PairType> >::
						 iterator ord_it = order.begin(); ord_it != order.end(); ++ord_it)
			{
				result.push_back(ord_it->second);
			}
		}

		/**
			 @brief Compute overall peptide abundances.

			 Based on quantitative data for individual charge states (derived from annotated features) in member @p pep_quant_, compute overall abundances for all peptides and store them also in @p pep_quant_.
		*/
		void quantifyPeptides_();

		/**
			 @brief Normalize peptide abundances across samples by (multiplicative) scaling to equal medians.
		*/
		void normalizePeptides_();

		/**
			 @brief Get the "canonical" protein accession from the list of protein accessions of a peptide.

			 @param pep_accessions Protein accessions of a peptide
			 @param accession_to_leader Captures information about indistinguishable proteins (maps accession to accession of group leader)

			 If there is no information about indistinguishable proteins (from protXML) available, a canonical accession exists only for proteotypic peptides - it's the single accession for the respective peptide.

			 Otherwise, a peptide has a canonical accession if it maps only to proteins of one indistinguishable group. In this case, the canonical accession is that of the group leader.

			 If there is no canonical accession, the empty string is returned.
		*/
		String getAccession_(const std::set<String>& pep_accessions, 
												 std::map<String, String>& accession_to_leader);

		/**
			 @brief Count the number of identifications (best hits only) of each peptide sequence.

			 The peptide hits in @p peptides are sorted by score in the process.
		*/
		void countPeptides_(std::vector<PeptideIdentification>& peptides);

		/// Clear all data when parameters are set
		void updateMembers_();

	}; // class

} // namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_PEPTIDEANDPROTEINQUANT_H
