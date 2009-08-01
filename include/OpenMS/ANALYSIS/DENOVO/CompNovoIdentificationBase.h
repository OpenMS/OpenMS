// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_DENOVO_COMPNOVOIDENTIFICATIONBASE_H
#define OPENMS_ANALYSIS_DENOVO_COMPNOVOIDENTIFICATIONBASE_H

// OpenMS includes
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
//#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include "MassDecomposition.h"
#include "MassDecompositionAlgorithm.h"
#include "CompNovoIonScoringBase.h"

// stl includes
#include <vector>

namespace OpenMS
{
	/**
	  @brief  run with CompNovoIdentificationBase

		@ref CompNovoIdentificationBase_Parameters are explained on a separate page.
		
		@ingroup Analysis_ID
	*/
	class CompNovoIdentificationBase : public DefaultParamHandler
	{

		public:

			/** @name constructors and destructors
			 */
			//@{
			/// default constructor
			CompNovoIdentificationBase();
			
			/// copy constructor
			CompNovoIdentificationBase(const CompNovoIdentificationBase& source);
			
			/// destructor
			virtual ~CompNovoIdentificationBase();
			//@}
		
			///
			CompNovoIdentificationBase& operator = (const CompNovoIdentificationBase& source);

			/** @name Accessors
			 */
			//@{
			/// performs an ProteinIdentification run on a PeakMap
			virtual void getIdentifications(std::vector<PeptideIdentification>& ids, const PeakMap& exp) = 0;

			/// performs an ProteinIdentification run on a PeakSpectrum
			//virtual void getIdentification(PeptideIdentification& id, const PeakSpectrum& CID_spec, const PeakSpectrum& ETD_spec) = 0;
			//@}

			typedef CompNovoIonScoringBase::IonScore IonScore;
			
		protected:

			/// initialize the needed members
			void init_();
			
			/// update members method from DefaultParamHandler to update the members 
			void updateMembers_();

			/// filters the permutations 
			void filterPermuts_(std::set<String>& permut);

			/// selects pivot ion of the given range using the scores given in CID_nodes
			void selectPivotIons_(std::vector<UInt>& pivots, UInt left, UInt right, Map<double, IonScore>& CID_nodes, const PeakSpectrum& CID_orig_spec, double precursor_weight, bool full_range = false);
		
			/// filters the decomps by the amino acid frequencies
			void filterDecomps_(std::vector<MassDecomposition>& decomps);
		
			/// produces mass decompositions using the given mass
			void getDecompositions_(std::vector<MassDecomposition>& decomps, double mass, bool no_caching = false);

			/// permuts the String s adds the prefix and stores the results in permutations
			void permute_(String prefix, String s, std::set<String>& permutations);
					
			UInt countMissedCleavagesTryptic_(const String& peptide) const;
			
			/// fills the spec with b and y ions, no other ion types or doubly charged variants are used
			void getCIDSpectrumLight_(PeakSpectrum& spec, const String& sequence, double prefix, double suffix);
			
			/// fills the spectrum with b,y ions, multiple charged variants; if prefix and suffix weights are given, the sequence is treated as tag
			void getCIDSpectrum_(PeakSpectrum& spec, const String& sequence, int charge, double prefix = 0.0, double suffix = 0.0);
		
			/// initializes the score distribution precalculated for the use in spectrum generation
			void initIsotopeDistributions_();

			/// estimates an exact precursor weight of the ETD spectrum, because in most of the cases the precursor is found in the MS/MS spec
			double estimatePrecursorWeight_(const PeakSpectrum& ETD_spec, UInt& charge);

			/// keep for each window of size windowsize in the m/z range of the spectrum exactly no_peaks
			void windowMower_(PeakSpectrum& spec, double windowsize, UInt no_peaks);

			/// compares two spectra 
			double compareSpectra_(const PeakSpectrum& s1, const PeakSpectrum& s2);

			/// returns a modified AASequence from a given internal representation
			AASequence getModifiedAASequence_(const String& sequence);

			/// returns the internal representation of a given AASequence
			String getModifiedStringFromAASequence_(const AASequence& sequence);

			/// mapping for the internal representation character to the actual residue
			Map<char, const Residue*> name_to_residue_;

			/// mapping of the actual residue to the internal representing character
			Map<const Residue*, char> residue_to_name_;
			
			///
			Map<int, std::vector<double> > isotope_distributions_;

			/// masses of the amino acids
			Map<char, double> aa_to_weight_; 

			MassDecompositionAlgorithm mass_decomp_algorithm_;

			double min_aa_weight_;

			ZhangSimilarityScore zhang_;

			Map<UInt, Map<UInt, std::set<String> > > subspec_to_sequences_;

			UInt max_number_aa_per_decomp_;

			bool tryptic_only_;

			double fragment_mass_tolerance_;

			UInt max_number_pivot_;

			double decomp_weights_precision_;

			double max_mz_;

			double min_mz_;

			double max_decomp_weight_;

			UInt max_subscore_number_;

			UInt max_isotope_;

			Map<double, std::vector<MassDecomposition> > decomp_cache_;

			Map<String, std::set<String> > permute_cache_;

			//PILISModel pilis_model_;
			
		public:
			
			class Permut
			{
				//public:
				private:
								
				Permut()
					: score_(0)
				{
				}

				public:
				
				Permut(const std::set<String>::const_iterator& permut, double s)
					: permut_(permut),
						score_(s)
				{
				}

				Permut(const Permut& rhs)
					: permut_(rhs.permut_),
						score_(rhs.score_)
				{	
				}

				~Permut()
				{
				}

				Permut& operator = (const Permut& rhs)
				{
					if (&rhs != this)
					{
						permut_ = rhs.permut_;
						score_ = rhs.score_;
					}
					return *this;
				}

				const std::set<String>::const_iterator& getPermut() const
				{
					return permut_;
				}

				void setPermut(const std::set<String>::const_iterator& it)
				{
					permut_ = it;
				}

				double getScore() const
				{
					return score_;
				}

				void setScore(double score)
				{
					score_ = score;
				}

			protected:

				std::set<String>::const_iterator permut_;
				double score_;

				//private:
					
				//Permut(); // private to prevent usage
			};

	};

	namespace Internal
	{
		bool PermutScoreComparator(const CompNovoIdentificationBase::Permut& p1, const CompNovoIdentificationBase::Permut& p2);
	}
}

#endif
