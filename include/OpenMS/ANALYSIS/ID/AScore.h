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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef  OPENMS_ANALYSIS_ID_ASCORE_H
#define  OPENMS_ANALYSIS_ID_ASCORE_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <vector>

namespace OpenMS
{
	class PeptideHit;
	class AASequence;
	struct ProbablePhosphoSites
	{
		Size first;
		Size second;
		Size seq_1;
		Size seq_2;
		Size peak_depth;
		Size AScore;
	};
	/**
		@brief Implementation of the Ascore
		For a given Peptidesequence and its MS^2 spectrum it is tried to identify the most probable phosphorylation-site(s). 
		For each phosphorylation site a score is calculated, which is an indicator for the probability that this site is phosphorylated. 
		The algorithm is implemented according to Beausoleil et al.	
		
	*/
	class OPENMS_DLLAPI AScore
	{
		public:
			///Default constructor
			AScore();
			///Destructor
			~AScore();

			/**
				@brief Computes the AScore and returns all computed phospho-sites. The safed sequences contain only phospho informations. All other modifications are dropped due to simplicity.

				@param	hit a PeptideHit
				@param real_spectrum spectrum mapped to hit
				@param fmt fragment_mass_tolarence, when comparing real_spectrum to a theoretical spectrum of the amino acid seqeunce of hit.
				@param number_of_phospho_sites which directs the method to search for this number of phosphorylated sites.
				
				@note the original sequence is safed in the PeptideHits as MetaValue Search_engine_sequence.
			*/
			PeptideHit compute(PeptideHit& hit, RichPeakSpectrum& real_spectrum, DoubleReal fmt, Int number_of_phospho_sites);
			
			///Computes the cumulative binomial probabilities.
			DoubleReal computeCumulativeScore(UInt N,UInt n, DoubleReal p);
			
			/**
				@brief Finds the peptides with the highest PeptideScores and outputs all informations for computing the AScore
				@note This function assumes that there are more permutations than the assumed number of phosphorylations!
			*/
			void computeHighestPeptides( std::vector< std::vector<DoubleReal> >& peptide_site_scores,std::vector<ProbablePhosphoSites>& sites, std::vector<std::vector<Size> >& permutations);
			///Computes the site determing_ions for the given AS and sequences in candidates
			void compute_site_determining_ions(std::vector<RichPeakSpectrum>& th_spectra, ProbablePhosphoSites& candidates, Int charge, std::vector<RichPeakSpectrum>& site_determining_ions);
		private:
			///computes number of matched ions between windows and the given spectrum. All spectra have to be sorted by Position!
			Int numberOfMatchedIons_(const RichPeakSpectrum& th,const RichPeakSpectrum& windows ,Size depth, DoubleReal fmt);
			///computes the peptide score according to Beausoleil et al. page 1291
			DoubleReal peptideScore_(std::vector<DoubleReal>& scores);
			public:
			///helperfunction 
			std::vector<Size> computeTupel_(AASequence& without_phospho);
			///helper function
			std::vector<std::vector<Size> > computePermutations_(std::vector<Size> tupel,Int number_of_phospho_sites);
	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_ASCORE_H
