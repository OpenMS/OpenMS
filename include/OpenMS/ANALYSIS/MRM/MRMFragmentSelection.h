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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MRM_MRMFRAGMENTSELECTION_H
#define OPENMS_ANALYSIS_MRM_MRMFRAGMENTSELECTION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief This class can select appropriate fragment ions of an MS/MS spectrum of a peptide

		@htmlinclude OpenMS_MRMFragmentSelection.parameters

		Several user choices can influence the selection of the ions from the MS/MS spectrum. These
		choices can be done using the parameters as described on the parameters page (see below). 
		Basically there are two different ways of selecting suitable ions. One, using standardized 
		names, e.g. given in the meta value "IonName" of each peaks of the spectrum (this can be 
		written from TheoreticalSpectrumGenerator, PILISModel...). The second one is simply using
		the most abundant peaks in a specified m/z range.

		@ingroup Analysis_MRM
	*/
	class OPENMS_DLLAPI MRMFragmentSelection : public DefaultParamHandler
	{

		public:

		/** @name Constructors and destructors
		*/
		//@{
		/// default constructor
		MRMFragmentSelection();

		/// copy constructor
		MRMFragmentSelection(const MRMFragmentSelection& rhs);

		/// destructor
		virtual ~MRMFragmentSelection();
		//@}

		/// assignment operator 
		MRMFragmentSelection& operator = (const MRMFragmentSelection& rhs);

		/// selects accordingly to the parameters the best peaks of spec and writes them into selected_peaks
		void selectFragments(std::vector<RichPeak1D>& selected_peaks, const RichPeakSpectrum& spec);
	
		protected:

		/// returns true if the selection of peak is allowed, according to the parameters set
		bool peakselectionIsAllowed_(const RichPeak1D& peak);
	};
}

#endif

