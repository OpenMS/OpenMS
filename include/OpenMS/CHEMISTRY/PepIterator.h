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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------
 
#ifndef OPENMS_CHEMISTRY_PEPITERATOR_H
#define OPENMS_CHEMISTRY_PEPITERATOR_H

#include <OpenMS/DATASTRUCTURES/String.h>
 
namespace OpenMS
{
 
	/**
	@brief Abstract base class for different peptide iterators
	
	@note every derived class has to implement the static functions "PepIterator * create()" and "const String getProductName()" (see Factory for details)
	*/
	class OPENMS_DLLAPI PepIterator
	{
 	
		public:
				
			typedef std::pair <String, String> FASTAEntry;
			
			/**
			@brief constructor
			*/
			PepIterator();
			
			/**
			@brief destructor
			*/
			virtual ~PepIterator();
			
			/**
			@brief copy constructor
			*/
			PepIterator(const PepIterator & source);
		
			/**
			@brief * operator for accessing the value of the iterator
			@return FASTAEntry representing a peptide
			@throw InvalidIterator if iterator has not been initialized
			*/
			virtual FASTAEntry operator*() = 0;
			
			/**
			@brief opperator ++ for postincrement
			@return Reference to PepIterator
			@throw InvalidIterator if iterator has not been initialized
			*/		
			virtual PepIterator & operator++() = 0;
			
			/**
			@brief opperator ++ for preincrement
			@return pointer to PepIterator
			@throw Exception::InvalidIterator if iterator has not been initialized
			*/
			virtual PepIterator * operator++(int) = 0;
		
			/**
			@brief setter for FASTA file
			@param f const String reference representing file location
			@throw Exception::FileNotFound
			@throw Exception::ParseError
			*/
			virtual void setFastaFile (const String & f) = 0;
			
			/**
			@brief getter for FASTA file
			@return String with file location
			*/
			virtual String getFastaFile() = 0;
			
			/**
			@brief setter for spectrum
			@param s ms spectrum given as vector of DoubleReals
			@throw Exception::InvalidValue if spectrum is not sorted acendingly
			*/
			virtual void setSpectrum (const std::vector<DoubleReal> & s) = 0;
		 
			/**
			@brief getter for spectrum
			@return the used spectrum
			*/
			virtual const std::vector<DoubleReal> & getSpectrum ()=0;
		
			/**
			@brief setter for tolerance
			@param t tolerance value
			@throw Exception::InvalidValue if tolerance is negative
			*/
			virtual void setTolerance (DoubleReal t) = 0;
			
			/**
			@brief getter for tolerance
			@return tolerance
			*/
			virtual DoubleReal getTolerance()=0;
			
			/**
			@brief initializing iterator
			*/
			virtual bool begin()=0;
			
			/**
			@brief idicator where iterator is at end
			*/
			virtual bool isAtEnd()=0;
			
			/**
				@brief all children has to be registered here
				@see Factory
			*/
			static void registerChildren();
		
	};
}
#endif // OPENMS_CHEMISTRY_PEPITERATOR_H
