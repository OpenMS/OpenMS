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


#ifndef OPENMS_ANALYSIS_DENOVO_COMPNOVOIONSCORINGBASE_H
#define OPENMS_ANALYSIS_DENOVO_COMPNOVOIONSCORINGBASE_H

// OpenMS includes
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

// stl includes
#include <vector>

namespace OpenMS
{
	/**
	  @brief  run with CompNovoIonScoringBase

		@ref CompNovoIonScoringBase_Parameters are explained on a separate page.
		
		@ingroup Analysis_ID
	*/
	class CompNovoIonScoringBase : public DefaultParamHandler
	{

		public:

			enum IsotopeType
			{
				PARENT = 0,
				CHILD = 1,
				LONE = 2
			};

			struct IonScore
			{
				IonScore();

				IonScore(const IonScore& rhs);

				virtual ~IonScore();

				IonScore& operator = (const IonScore& rhs);

				
				double score;
        double s_bion;
        double s_yion;
        double s_witness;
        double position;
        double s_isotope_pattern_1; // isotope pattern score charge 1
        int is_isotope_1_mono; // 0 means not testet, 1 mean is, -1 is tail of isotopes
        double s_isotope_pattern_2; // "" charge 2
			};

						
			/** @name constructors and destructors
			 */
			//@{
			/// default constructor
			CompNovoIonScoringBase();
			
			/// copy constructor
			CompNovoIonScoringBase(const CompNovoIonScoringBase& source);
			
			/// destructor
			virtual ~CompNovoIonScoringBase();
			//@}
		
			///
			CompNovoIonScoringBase& operator = (const CompNovoIonScoringBase& source);

			/** @name Accessors
			 */
			//@{
			double scoreIsotopes(const PeakSpectrum& CID_spec, PeakSpectrum::ConstIterator it, UInt charge);
			//@}

		protected:

			/// update members method from DefaultParamHandler to update the members 
			void updateMembers_();


			IsotopeType classifyIsotopes_(const PeakSpectrum& spec, PeakSpectrum::ConstIterator it);

			double scoreIsotopes_(const PeakSpectrum& spec, PeakSpectrum::ConstIterator it, Map<double, IonScore>& CID_nodes, UInt charge = 1);

			virtual void scoreWitnessSet_(UInt charge, double precursor_weight, Map<double, IonScore>& CID_nodes, const PeakSpectrum& CID_orig_spec) = 0;
			
			void addSingleChargedIons_(Map<double, IonScore>& ion_scores, PeakSpectrum& CID_spec);

			void initIsotopeDistributions_();

			///
			Map<int, std::vector<double> > isotope_distributions_;

			double fragment_mass_tolerance_;
		
		public:

	};

}

#endif
