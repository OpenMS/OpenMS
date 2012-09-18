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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_ID_PILISIDENTIFICATION_H
#define OPENMS_ANALYSIS_ID_PILISIDENTIFICATION_H

#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>
#include <map>

namespace OpenMS
{
	// forward declarations
	class PeakSpectrumCompareFunctor;
	
	/**
	  @brief This class actually implements a complete ProteinIdentification run with PILIS

		The PILISIdentification class needs a PILISModel and a PILISSequenceDB to generate
		identifications. Simply call getIdentifications with a RichPeakMap.
		 
		@htmlinclude OpenMS_PILISIdentification.parameters
		
		@ingroup Analysis_ID
	*/
	class OPENMS_DLLAPI PILISIdentification : public DefaultParamHandler
	{

		public:

			/** @name constructors and destructors
			 */
			//@{
			/// default constructor
			PILISIdentification();
			
			/// copy constructor
			PILISIdentification(const PILISIdentification& source);
			
			/// destructor
			virtual ~PILISIdentification();
			//@}
		
			///
			PILISIdentification& operator = (const PILISIdentification& source);

			/** @name Accessors
			 */
			//@{
			/// sets the sequence DB to be used for the ProteinIdentification runs
			//void setSequenceDB(const SuffixArrayPeptideFinder& sapf);

			/// sets the model to be used for the ProteinIdentification run
			void setModel(PILISModel* hmm_model);

			/// performs an ProteinIdentification run on a RichPeakMap
			void getIdentifications(const std::vector<std::map<String, UInt> >& candidates, std::vector<PeptideIdentification>& ids, const RichPeakMap& exp);

			/// performs an ProteinIdentification run on a PeakSpectrum
			void getIdentification(const std::map<String, UInt>& candidates, PeptideIdentification& id, const RichPeakSpectrum& spectrum);
			//@}

		protected:

			/// fast method to create spectra for pre-scoring
			void getSpectrum_(RichPeakSpectrum& spec, const String& sequence, int charge);
		
			/// performs a pre-scoring of the given spec with very simple spectra from the candidate peptides
			void getPreIdentification_(PeptideIdentification& id, const RichPeakSpectrum& spec, const std::map<String, UInt>& cand_peptides);

			/// performs a ProteinIdentification via spectra comparison with the PILISModel spectrum generator
			void getFinalIdentification_(PeptideIdentification& id, const RichPeakSpectrum& spec, const PeptideIdentification& pre_id);
	
			/// returns the model pointer
			PILISModel* getPILISModel_();

			/// returns the sequence database pointer
			//PILISSequenceDB* getSequenceDB_();
			//SuffixArrayPeptideFinder* getSequenceDB_();
			
			/// the sequence database for the candidate peptides
			//PILISSequenceDB* sequence_db_;
			//SuffixArrayPeptideFinder sapf_;

			/// the model for spectra simulation
			PILISModel* hmm_model_;

			/// amino acids weights for the simple spectra generator
			Map<char, double> aa_weight_;

			/// scorer for pre comparison
			PeakSpectrumCompareFunctor* pre_scorer_;

			/// scorer for spectra comparison
			PeakSpectrumCompareFunctor* scorer_;

			/// a peaks, just to not instantiate it over and over again
			RichPeak1D p_;

			///
			std::vector<RichPeakSpectrum> sim_specs_;

			/// flag whether the istance has a internal sequence db
			bool own_sequence_db_;

			/// flag whether the istance has a internal model
			bool own_model_;

			/// update members method from DefaultParamHandler to update the members 
			void updateMembers_();
	};
}

#endif
