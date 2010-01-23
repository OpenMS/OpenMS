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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATORTRYPTIC_H
#define OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATORTRYPTIC_H

#include <OpenMS/CHEMISTRY/EdwardsLippertIterator.h>

namespace OpenMS 
{
/**
@brief EdwardsLippertIterator that only retrieves tryptic seqences
*/
class OPENMS_DLLAPI EdwardsLippertIteratorTryptic : public EdwardsLippertIterator
{
	public:

		/// default constructor
		EdwardsLippertIteratorTryptic();

		/// copy constructor
		EdwardsLippertIteratorTryptic(const EdwardsLippertIteratorTryptic& rhs);

		/// destructor
		virtual ~EdwardsLippertIteratorTryptic();
					
		/// assignment operator
		EdwardsLippertIteratorTryptic& operator = (const EdwardsLippertIteratorTryptic& rhs);
		
		/**
		@brief indicates if trypsin will cat between the two amino acids
		@param aa1 first amino acid
		@param aa2 second amino acid
		*/
		virtual bool isDigestingEnd (char aa1, char aa2);
	
		/**
		@brief needed by Factory
		@return const string name of class
		*/
		static const String getProductName()
		{
			return "EdwardsLippertIteratorTryptic";
		}

		/**
		@brief needed by Factory
		@return poiter to new object
		*/
		static PepIterator * create()
		{
			return new EdwardsLippertIteratorTryptic;
		}
	};

}

#endif //OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATORTRYPTIC_H
