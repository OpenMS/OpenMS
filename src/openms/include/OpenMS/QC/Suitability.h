// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <map>
#include <vector>

namespace OpenMS
{
	class MSExperiment;
	class PeptideIdentification;
	class PeptideHit;

	/**
	* @brief This class serves as the library representation of @ref TOPP_DatabaseSuitability
	*
	* This class holds the functionality of calculating the database suitability and everything
	* that is needed for this. Quality of LC-MS/MS spectra can also be calculated.
	*/
	class OPENMS_DLLAPI Suitability
	{
	public:
		/// Constructor
		Suitability() = default;

		/// Destructor
		~Suitability() = default;

		/**
		* @brief Computes quality of LC-MS/MS spectra (id-rate of deNovo sequences)
		*
		* spectra quality = #deNovo seqs / #MS2 spectra
		*
		* @param exp MSExperiment from which the deNovo sequences where calculated
		* @param pep_ids vector containing the deNovo sequences as pepIDs
		* @return spectral quality
		* @throws MissingInformation if no MS2 spectra are found
		*/
		double computeSpectraQuality(const MSExperiment& exp, const std::vector<PeptideIdentification>& pep_ids);

		double computeSuitability(const std::vector<PeptideIdentification>& pepIDs, double FDR, double novo_fract, bool no_re_rank);

		std::map<String, double> getData();

	private:
		Size num_novo_seqs;
		Size num_ms2;
		Size num_unique_novo_seqs;
		double spectral_quality;
		Size num_top_novo;
		Size num_top_db;
		Size num_re_ranked;
		Size num_interest;
		double cut_off;
		double suitability;

		/**
		* @brief Calculates the xcorr difference between the top two hits marked as decoy
		*
		* Only searches the top ten hits for two decoys. If there aren't two decoys, DBL_MAX
		* is returned.
		*
		* @param pep_id pepID from where the decoy difference will be calculated
		* @return xcorr difference
		* @throws MissingInformation if no target/decoy annotation is found
		* @throws MissingInformation if no xcorr is found
		*/
		double getDecoyDiff_(const PeptideIdentification& pep_id);

		/**
		* @brief Calculates a xcorr cut-off based on decoy hits
		*
		* All N decoy differences are calculated. The (1-novor_fract)*N highest one
		* is returned.
		* It is asssumed that this difference accounts for novor_fract re-ranking cases.
		*
		* @param pep_ids vector containing the pepIDs
		* @param novor_fract fraction of re-ranking cases should be re-ranked
		* @return xcorr cut-off
		* @throws MissingInformation if no more than 20 % of the pepIDs have two decoys in there top ten
		*/
		double getDecoyCutOff_(const std::vector<PeptideIdentification>& pep_ids, double novor_fract);

		/**
		* @brief Tests if a PeptideHit is considered a deNovo hit
		*
		* To test this the function looks into the protein accessions.
		* If only the deNovo protein is found, 'true' is returned.
		* If one database protein is found, 'false' is returned.
		*
		* @param hit PepHit in question
		* @return true/false
		*/
		bool isNovoHit_(const PeptideHit& hit);

		/**
		* @brief Tests if a PeptideHit has a higher q-value the given FDR
		*
		* Q-value is searched at score and at meta-value level.
		*
		* @param hit PepHit in question
		* @param FDR value to check against
		* @param q_value_score is q-value the current score type?
		* @return true/false
		* @throws Precondition if no q-value is found
		*/
		bool scoreHigherThanFDR_(const PeptideHit& hit, double FDR, bool q_value_score);
	};
}

