// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_DENOVO_COMPNOVOIDENTIFICATIONBASE_H
#define OPENMS_ANALYSIS_DENOVO_COMPNOVOIDENTIFICATIONBASE_H

// OpenMS includes
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringBase.h>

// stl includes
#include <vector>

namespace OpenMS
{
	/**
	  @brief  run with CompNovoIdentificationBase

		@ingroup Analysis_ID
	*/
	class OPENMS_DLLAPI CompNovoIdentificationBase : public DefaultParamHandler
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
		
			/// assignment operator
			CompNovoIdentificationBase& operator = (const CompNovoIdentificationBase& source);

			/** @name Accessors
			 */
			//@{
			/// performs an ProteinIdentification run on a PeakMap
			virtual void getIdentifications(std::vector<PeptideIdentification>& ids, const PeakMap& exp) = 0;
			//@}

			typedef CompNovoIonScoringBase::IonScore IonScore;
			
		protected:

			/// update members method from DefaultParamHandler to update the members 
			void updateMembers_();

			/// filters the permutations 
			void filterPermuts_(std::set<String>& permut);

			/// selects pivot ion of the given range using the scores given in CID_nodes
			void selectPivotIons_(std::vector<Size>& pivots, Size left, Size right, Map<DoubleReal, IonScore>& CID_nodes, const PeakSpectrum& CID_orig_spec, DoubleReal precursor_weight, bool full_range = false);
		
			/// filters the decomps by the amino acid frequencies
			void filterDecomps_(std::vector<MassDecomposition>& decomps);
		
			/// produces mass decompositions using the given mass
			void getDecompositions_(std::vector<MassDecomposition>& decomps, DoubleReal mass, bool no_caching = false);

			/// permuts the String s adds the prefix and stores the results in permutations
			void permute_(String prefix, String s, std::set<String>& permutations);
					
			Size countMissedCleavagesTryptic_(const String& peptide) const;
			
			/// fills the spec with b and y ions, no other ion types or doubly charged variants are used
			void getCIDSpectrumLight_(PeakSpectrum& spec, const String& sequence, DoubleReal prefix, DoubleReal suffix);
			
			/// fills the spectrum with b,y ions, multiple charged variants; if prefix and suffix weights are given, the sequence is treated as tag
			void getCIDSpectrum_(PeakSpectrum& spec, const String& sequence, Size charge, DoubleReal prefix = 0.0, DoubleReal suffix = 0.0);
		
			/// initializes the score distribution precalculated for the use in spectrum generation
			void initIsotopeDistributions_();

			/// estimates an exact precursor weight of the ETD spectrum, because in most of the cases the precursor is found in the MS/MS spec
			DoubleReal estimatePrecursorWeight_(const PeakSpectrum& ETD_spec, Size& charge);

			/// keep for each window of size windowsize in the m/z range of the spectrum exactly no_peaks
			void windowMower_(PeakSpectrum& spec, DoubleReal windowsize, Size no_peaks);

			/// compares two spectra 
			DoubleReal compareSpectra_(const PeakSpectrum& s1, const PeakSpectrum& s2);

			/// returns a modified AASequence from a given internal representation
			AASequence getModifiedAASequence_(const String& sequence);

			/// returns the internal representation of a given AASequence
			String getModifiedStringFromAASequence_(const AASequence& sequence);

			/// mapping for the internal representation character to the actual residue
			Map<char, const Residue*> name_to_residue_;

			/// mapping of the actual residue to the internal representing character
			Map<const Residue*, char> residue_to_name_;
			
			///
			Map<Size, std::vector<DoubleReal> > isotope_distributions_;

			/// masses of the amino acids
			Map<char, DoubleReal> aa_to_weight_; 

			MassDecompositionAlgorithm mass_decomp_algorithm_;

			DoubleReal min_aa_weight_;

			ZhangSimilarityScore zhang_;

			Map<Size, Map<Size, std::set<String> > > subspec_to_sequences_;

			Size max_number_aa_per_decomp_;

			bool tryptic_only_;

			DoubleReal fragment_mass_tolerance_;

			Size max_number_pivot_;

			DoubleReal decomp_weights_precision_;

			DoubleReal max_mz_;

			DoubleReal min_mz_;

			DoubleReal max_decomp_weight_;

			Size max_subscore_number_;

			Size max_isotope_;

			Map<DoubleReal, std::vector<MassDecomposition> > decomp_cache_;

			Map<String, std::set<String> > permute_cache_;

		public:
			
			/** @brief Simple class to store permutations and a score

					This class is used to store the generated perumtations
					and a score to them
			*/
			class Permut
			{
				private:
								
				Permut()
					: score_(0)
				{
				}

				public:
				
				Permut(const std::set<String>::const_iterator& permut, DoubleReal s)
					: permut_(permut),
						score_(s)
				{
				}

				Permut(const Permut& rhs)
					: permut_(rhs.permut_),
						score_(rhs.score_)
				{	
				}

				virtual ~Permut()
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

				DoubleReal getScore() const
				{
					return score_;
				}

				void setScore(DoubleReal score)
				{
					score_ = score;
				}

			protected:

				std::set<String>::const_iterator permut_;
				DoubleReal score_;
			};

	};

	namespace Internal
	{
		bool PermutScoreComparator(const CompNovoIdentificationBase::Permut& p1, const CompNovoIdentificationBase::Permut& p2);
	}
}

#endif
