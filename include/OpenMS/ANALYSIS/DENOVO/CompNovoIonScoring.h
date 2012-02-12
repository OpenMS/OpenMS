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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_DENOVO_COMPNOVOIONSCORING_H
#define OPENMS_ANALYSIS_DENOVO_COMPNOVOIONSCORING_H

// OpenMS includes
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
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
	  @brief  run with CompNovoIonScoring

		@htmlinclude OpenMS_CompNovoIonScoring.parameters
		
		@ingroup Analysis_DeNovo
	*/
	class OPENMS_DLLAPI CompNovoIonScoring : public CompNovoIonScoringBase
	{

		public:

			typedef CompNovoIonScoringBase::IsotopeType IsotopeType;
			typedef CompNovoIonScoringBase::IonScore IonScore;

						
			/** @name constructors and destructors
			 */
			//@{
			/// default constructor
			CompNovoIonScoring();
			
			/// copy constructor
			CompNovoIonScoring(const CompNovoIonScoring& source);
			
			/// destructor
			virtual ~CompNovoIonScoring();
			//@}
		
			/// assignment operator
			CompNovoIonScoring& operator = (const CompNovoIonScoring& source);

			/** @name Accessors
			 */
			//@{
			void scoreSpectra(Map<DoubleReal, IonScore>& CID_ion_scores, PeakSpectrum& CID_spec, PeakSpectrum& ETD_spec, DoubleReal precursor_weight, Size charge);
			//@}

		protected:

			void scoreETDFeatures_(Size charge, DoubleReal precursor_weight, Map<DoubleReal, IonScore>& CID_nodes, const PeakSpectrum& CID_orig_spec, const PeakSpectrum& ETD_orig_spec);

			void scoreWitnessSet_(Size charge, DoubleReal precursor_weight, Map<DoubleReal, IonScore>& CID_nodes, const PeakSpectrum& CID_orig_spec);
			
	};

}

#endif
