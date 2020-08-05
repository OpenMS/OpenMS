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

#include <vector>

namespace OpenMS
{
	class MSExperiment;
	class PeptideIdentification;

	/**
	* @brief Spectral quality based on deNovo identification rate
	*
	* Simple class the calculate the id-rate of a vector of pepIDs
	* This can be used to calculate spectral quality when the pepIDs
	* are deNovo sequences calculated from the mzML.
	* 
	* TODO: Change MS2IdentificationRate so it can be used instead
	*
	*/
	class OPENMS_DLLAPI SpectralQuality
	{
	public:

		struct SpectralData
		{
			Size num_novo_seqs = 0;
			Size num_ms2 = 0;
			Size num_unique_novo_seqs = 0;
			double spectral_quality = 0;
		};

		SpectralQuality() = default;

		~SpectralQuality() = default;

		/**
			* @brief Computes quality of LC-MS/MS spectra (id-rate of deNovo sequences)
			*
			* spectra quality = #deNovo seqs / #MS2 spectra
			* The results are written into a SpectraData object which is then appended to the result vector
			*
			* @param exp				MSExperiment from which the deNovo sequences where calculated
			* @param pep_ids		vector containing the deNovo sequences as pepIDs
			* @throws						MissingInformation if no MS2 spectra are found
			*/
		void computeSpectraQuality(const MSExperiment& exp, const std::vector<PeptideIdentification>& pep_ids);

		std::vector<SpectralData> getResults() const;

	private:
		std::vector<SpectralData> results;
	};
}

