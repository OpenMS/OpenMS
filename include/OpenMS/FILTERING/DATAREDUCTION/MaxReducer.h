// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_MAXREDUCER_H
#define OPENMS_FILTERING_DATAREDUCTION_MAXREDUCER_H

#include<OpenMS/FILTERING/DATAREDUCTION/DataReducer.h>

namespace OpenMS
{	
	/**
		@brief Reduces the amount of data in an experiment by extracting the maximum intensity peaks.
		
	*/
	class MaxReducer
	  : public DataReducer
	{
		public:

			///constructor
			MaxReducer();
	
			///destructor
			~MaxReducer();
	
			/// calculates the reduced MSExperiment
		 	virtual  void applyReduction(const  ExperimentType& in,  ExperimentType& out );
		
			/// returns an instance of this class
			static DataReducer* create()
			{
				return new MaxReducer();
			}
	
			/// returns the name of this module
			static const String getProductName()
			{
				return "MaxReducer";
			}
	};
}
#endif
