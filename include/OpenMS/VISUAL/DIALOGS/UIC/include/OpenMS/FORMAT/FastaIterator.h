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


#ifndef OPENMS_FORMAT_FASTAITERATOR_H
#define OPENMS_FORMAT_FASTAITERATOR_H

#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <fstream>
#include <vector>

namespace OpenMS
{

/**
@brief Iterator over FASTA file

iterates over FASTA file without loading it into memory. It just holds just one entry in memory.
@see FastaIteratorIntern.h
*/ 
class OPENMS_DLLAPI FastaIterator : public PepIterator
{
	
	public:
		
	typedef std::pair <String, String> FASTAEntry;
	/**
	@brief constructor
	*/
	FastaIterator();
	
	/**
	@brief copy constructor
	*/
	FastaIterator(const FastaIterator &);
	
	/**
	@brief destructor
	*/
	virtual ~FastaIterator ();

	/**
	@brief * operator for getting the iterator's value
	@return FASTAEntry
	@throw InvalidIterator if iterator was not initialized
	*/
	virtual FASTAEntry operator*();
	
	/**
	@brief postincrement Operator for the iterator 
	@return reference to PepIterator	
	@throw Exception::InvalidIterator if iterator was not initialized
	*/
	virtual PepIterator & operator++();
	
	/**
	@brief preincrement Operator for the iterator
	@return pointer to PepIterator
	@throw Exception::InvalidIterator if iterator was not initialized
	*/
	virtual PepIterator * operator++(int i);
	
	/**	
	@brief setter for FASTAfile
	@param f Name of the fasta file
	@throw Exception::FileNotFound 
	@throw ParseError is thrown if the file could not be parsed
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
	virtual void setSpectrum (const std::vector<DoubleReal> &)
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
	static const String getProductName()
	{
		return "FastaIterator";
	}
	
	/**
	@brief needed by Factory
	@return poiter to new object
	*/
	static PepIterator * create()
	{
		return new FastaIterator;
	}
	
	
	
	protected:
		/**
		@brief gets the next string
		@return string
		*/
		virtual std::string next_ ();

		bool is_at_end_; ///< bool indicated whether iterator is at end

		std::ifstream * input_file_; ///< input file

		String fasta_file_; ///< fasta file location

		std::string actual_seq_; ///< actual sequence

		std::string header_; ///< actual fasta header

		std::string last_header_; ///< last fasta header
	
};
}
#endif // OPENMS_FORMAT_FASTAITERATOR_H
