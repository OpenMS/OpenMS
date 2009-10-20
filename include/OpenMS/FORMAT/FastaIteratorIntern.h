// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FASTAITERATORINTERN_H
#define OPENMS_FORMAT_FASTAITERATORINTERN_H

#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <vector>

namespace OpenMS
{

	/**
		@brief Iterator for a FASTA file
		
		In comparision to FastaIterator the FASTA file will be loaded first and stored to RAM, while the FastaIterator just iterates over the FASTA file without loading it completly to memory.
		
		@see FastaIterator
	*/
	class OPENMS_DLLAPI FastaIteratorIntern : public PepIterator
	{
	
		public: 
			
			typedef std::pair<String,String> FASTAEntry;
			
			/**
			@brief constructor
			*/
			FastaIteratorIntern();
			
			/**
			@brief copy constructor
			*/
			FastaIteratorIntern(const FastaIteratorIntern&);
		
			/**
			@brief constructor
			*/
			virtual ~FastaIteratorIntern();
		
			/**
			@brief * Operator for derefering of iterator
			@return FASTEEntry iterator at actual position
			@throw Exception::InvalidIterator if iterator has not been initialized
			*/
			virtual FASTAEntry operator*();
			
			/**
			@brief ++ Operator for the iterator
			@return reference to PepIterator
			@throw Exception::InvalidIterator if iterator has not been initialized
			*/
			virtual PepIterator & operator++();
			
			/**
			@brief ++ Operator for the iterator
			@return pointer to PepIterator
			@throw Exception::InvalidIterator if iterator has not been initialized
			*/
			virtual PepIterator * operator++(int i);
			
			/**
			@brief setter for FASTA file
			@param f const String reference representing file location
			@throw Exception::FileNotFound
			@throw Exception::ParseError
			*/
			virtual void setFastaFile (const String & f);
		
			/**
			@brief getter for FASTA file
			@return String with file location
			*/
			virtual String getFastaFile();
			
			/**
			@brief setter for spectrum
			@note note availeble for FastaIterator
			@throw Exception::NotImplemented 
			*/
			virtual void setSpectrum (const std::vector<DoubleReal> & /*spec*/)
			{
				throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			};
		
			/**
			@brief getter for spectrum
			@note note availeble for FastaIterator
			@throw Exception::NotImplemented 
			*/
			virtual const std::vector<DoubleReal> & getSpectrum ()
			{
				throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			};
			
			/**
			@brief setter for tolerance
			@note note availeble for FastaIterator
			@throw Exception::NotImplemented 
			*/
			virtual void setTolerance (DoubleReal /* t */)
			{
				throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			};
			
			/**
			@brief getter for tolerance
			@note note availeble for FastaIterator
			@return tolerance
			@throw Exception::NotImplemented 
			*/
			virtual DoubleReal getTolerance()
			{
				throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			};
		
			/**
			@brief initializing of iterator
			@return true if everything went rigth
			@throw Exception::InvalidIterator if fastaFile was not set
			*/	
			virtual bool begin ();
			
			/**
			@brief indicates whether iterator is at end
			@return bool true if interator is at end
			*/
			virtual bool isAtEnd ();
		
			/**
			@brief needed by Factory
			@return const string name of class
			*/
			static const std::string getProductName()
			{
				return "FastaIteratorIntern";
			}
		
			/**
			@brief needed by Factory
			@return poiter to new object
			*/
			static PepIterator * create()
			{
				return new FastaIteratorIntern;
			}
			
		
		
		protected:
	
			String fasta_file_; ///< location of the fasta file
		
			std::vector<FASTAEntry > entrys_; ///< content of fasta file
		
			std::vector<FASTAEntry >::iterator it_; ///< iterator over fasta file content
	
	};

}

#endif //OpenMS/FORMAT/FastaIteratorIntern
