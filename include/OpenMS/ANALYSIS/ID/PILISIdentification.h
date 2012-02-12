// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
