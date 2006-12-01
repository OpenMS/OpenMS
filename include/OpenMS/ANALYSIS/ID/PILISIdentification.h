// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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


#ifndef OPENMS_ANALYSIS_ID_PILISIDENTIFICATION_H
#define OPENMS_ANALYSIS_ID_PILISIDENTIFICATION_H

#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/ANALYSIS/ID/PILISSequenceDB.h>

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

namespace OpenMS
{
	/**
	  @brief This class actually implements a complete identification run with PILIS

		The PILISIdentification class needs a PILISModel and a PILISSequenceDB to generate
		identifications. Simply call getIdentifications with a PeakMap.
	*/
	class PILISIdentification
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
			const PILISIdentification& operator = (const PILISIdentification&);

			/** @name Accessors
			 */
			//@{
			///
			void setSequenceDB(PILISSequenceDB* sequence_db);

			///
			void setModel(PILISModel* hmm_model);

			///
			void getIdentifications(std::vector<Identification>& ids, const PeakMap& exp);

			/// returns the parameters
			Param& getParam();

			/// sets the parameters
			void setParam(const Param& param);
			//@}

		protected:

			Param param_;

			PILISSequenceDB* sequence_db_;

			PILISModel* hmm_model_;
	};
}

#endif
