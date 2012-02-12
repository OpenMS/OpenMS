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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_GRADIENT_H
#define OPENMS_METADATA_GRADIENT_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Representation of a HPLC gradient
    
    It consists of several eluents and timepoints.
    Linear behaviour between timepoints is assumed.
		
		@ingroup Metadata
  */
  class OPENMS_DLLAPI Gradient
  {
		public:
      /// Constructor
      Gradient();
      /// Copy constructor
      Gradient(const Gradient& source);
      /// Destructor
      ~Gradient();
      
      /// Assignment operator
      Gradient& operator = (const Gradient& source);
      
      /// Equality operator
      bool operator == (const Gradient& source) const;
      /// Equality operator
      bool operator != (const Gradient& source) const;
			
			/**
				@brief Adds an eluent at the end of the eluent array
			
				@exception Exception::InvalidValue is thrown if the same eluent name is used twice.
			*/
			void addEluent(const String& eluent);
			/// removes all eluents
			void clearEluents();
			/// returns a const reference to the list of eluents
			const std::vector <String>& getEluents() const;
	
			/**
				@brief Adds a timepoint at the end of the timepoint array
			
				@exception Exception::OutOfRange is thrown if the new timpoint is before the last timepoint.
			*/
			void addTimepoint(Int timepoint);
			/// removes all timepoints
			void clearTimepoints();
			/// returns a const reference to the list of timepoints
			const std::vector<Int>& getTimepoints() const;		

			/**
				@brief sets the percentage of eluent @p eluent at timepoint @p timepoint
			
				@exception Exception::InvalidValue is thrown if the eluent, timepoint or percentage is invalid.
			*/
			void setPercentage(const String& eluent, Int timepoint, UInt percentage);
			
			/**
				@brief returns a const reference to the percentages
				
				First dimension of the vector is the eluents, second dimension is the timepoints.
			*/
			const std::vector< std::vector< UInt > >& getPercentages() const;
			
			/**
				@brief returns the percentage of an @p eluent at a @p timepoint
				
				@exception Exception::InvalidValue is thrown if the eluent or timepoint is invalid.
			*/
			UInt getPercentage(const String& eluent, Int timepoint) const;
			
			/// sets all precentage values to 0
			void clearPercentages();
			
			/// checks if the percentages of all timepoints add up to 100%
			bool isValid() const;
			
    protected:
      std::vector<String> eluents_;
      std::vector<Int> times_;
      // first dimension is eluents, second is times
      std::vector< std::vector< UInt > > percentages_;
  };
 
} // namespace OpenMS

#endif // OPENMS_METADATA_GRADIENT_H

