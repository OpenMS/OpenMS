// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_ANALYSIS_ID_PILISMODELGENERATOR_H
#define OPENMS_ANALYSIS_ID_PILISMODELGENERATOR_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS 
{
	class HiddenMarkovModel;

	/** 
	  @brief This class implements the simulation of the spectra from PILIS

		PILIS uses a HMM based structure to model the population of fragment ions
		from a peptide. The spectrum generator can be accessed via the getSpectrum
		method.

		@htmlinclude OpenMS_PILISModelGenerator.parameters
	*/	
	class OPENMS_DLLAPI PILISModelGenerator : public DefaultParamHandler
	{
	
		public:

			/** @name Constructors and destructors
			*/
			//@{
			/// default constructor
			PILISModelGenerator();

			/// copy constructor
			PILISModelGenerator(const PILISModelGenerator& model);
			
			/// destructor
			virtual ~PILISModelGenerator();
			//@}

			/// assignment operator
			PILISModelGenerator& operator = (const PILISModelGenerator& mode);
		

			/** @name Accessors
			*/
			//@{
			/// generates the model and writes it into model
			void getModel(HiddenMarkovModel& model);

			/// generates the precursor model and writes it into model
			//void getPrecursorModel(HiddenMarkovModel& model);
			//@}


		protected:

			/// initializes the model
			//void initModels_();

			///
			//void initMainModel_();

			///
			//void initLossModels_();

			///
			//void initPrecursorModel_();
	};
}
#endif

