// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Sandro Andreotti $
// $Authors: Sandro Andreotti, Andreas Bertsch$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DENOVO_DENOVOIONSCORING_H
#define OPENMS_ANALYSIS_DENOVO_DENOVOIONSCORING_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <vector>

namespace OpenMS
{
  /**
    @brief Base class for ion scoring implementation for de novo algorithms


    
		@ingroup Analysis_DeNovo
  */
  class OPENMS_DLLAPI DeNovoIonScoring
  	: public DefaultParamHandler
  {
  	public:	
      /** @brief IonScore

					IonScore describes the likelihood of position to be prefix residue masses (in case of 
					an b-ion centric algorithm). 
			*/
			class IonScore
			{
			 public:
        /// score 
        DoubleReal score;

        /// position of the ion
        DoubleReal position;

        /// index of peak in the spectrum, -1 if not in spectrum
        SignedSize index;


				IonScore()
					: score(0),
						position(0.0),
						index(-1)
				{
				}

				IonScore(const IonScore& rhs)
					: score(rhs.score),
						position(rhs.score),
						index(rhs.index)
				{
				}

				virtual ~IonScore()
				{
				}

				IonScore& operator = (const IonScore& rhs)
				{
					if (this != &rhs)
					{
						score = rhs.score;
						position = rhs.position;
						index = rhs.index;
					}
					return *this;
				}

			};


			/** @name Constructors and destructors
			*/
			//@{
	  	/// default constructor
	  	DeNovoIonScoring();
  	
			/// destructor
			virtual ~DeNovoIonScoring();

  		/// copy constructor
  		DeNovoIonScoring(const DeNovoIonScoring& rhs);
			//@}
  		
			/// assignment operator
			DeNovoIonScoring& operator = (const DeNovoIonScoring& rhs);
		
			virtual void getIonScores(std::vector<IonScore>& ion_scores, const RichPeakSpectrum& spec) = 0;

			virtual void getIonScores(std::vector<std::vector<IonScore> >& ion_scores, const RichPeakMap& exp) = 0; 
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_DENOVO_DENOVOIONSCORING_H
