// -*- mode: C++; tab-width: 2; -*-
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


#ifndef OPENMS_ANALYSIS_ID_PILISSCORING_H
#define OPENMS_ANALYSIS_ID_PILISSCORING_H

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

namespace OpenMS
{
	/**
	  @brief This class actually implements the E-value based scoring of PILIS

		The method which is used to generate the E-values are adapted from
	
		David Fenyo and Ronald C. Beavis
		Anal. Chem. 2003, 75, 768-774
		A Method for Assessing the Statistical Significance of Mass Spectrometry-Based Protein
		Identifications Using General Scoring Schemes.

		The bases for the calculation are the similarity scores of the simulated
		and experimental spectra. The scores are tranformed into a discrete 
		score distribution and from this distribution E-values are calculated for
		the peptide hits. 

		If more than one spectrum is given two E-values can be calculated, one which gives
		the significance of the peptide hit considering only one spectrum, and the other also
		considering also all other hits of all other spectra. The second type of scoring
		is somewhat more accurate.
		 
		@htmlinclude OpenMS_PILISScoring.parameters

		@ingroup Analysis_ID
	*/
	class OPENMS_DLLAPI PILISScoring : public DefaultParamHandler
	{

		public:

			/** @name constructors and destructors
			 */
			//@{
			/// default constructor
			PILISScoring();
			
			/// copy constructor
			PILISScoring(const PILISScoring& source);
			
			/// destructor
			virtual ~PILISScoring();
			//@}
		
			///
			PILISScoring& operator = (const PILISScoring& source);

			/** @name Accessors
			 */
			//@{
			/// performs an ProteinIdentification run on a PeakMap
			void getScores(std::vector<PeptideIdentification>& ids);

			/// performs an ProteinIdentification run on a PeakSpectrum
			void getScore(PeptideIdentification& id);
			//@}

		protected:

			///
			void getFitParameter_(double& slope, double& intercept, const std::vector<double>& scores, double threshold);

			///
			void getSurvivalFunction_(Map<UInt, double>& points, std::vector<DPosition<2> >& survival_function);

			///
			void getScore_(PeptideIdentification& id, double global_slope, double global_intercept);
	
	};
}

#endif
