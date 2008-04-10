// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bauer$
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_TRYPTICITERATOR_H
#define OPENMS_CHEMISTRY_TRYPTICITERATOR_H

#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>

namespace OpenMS {
/**
@brief finds all tryptic Peptides with every missed cleavage

*/
class TrypticIterator : public PepIterator
{
	
	public:
	
	typedef std::pair<String,String> FASTAEntry;	
	/**
	@brief Constructor
	*/
	TrypticIterator();
	/**
	@brief Copy Constructor
	*/
	TrypticIterator (const TrypticIterator &);
	/**
	@brief Destructor
	*/
	virtual ~TrypticIterator ();
	
	/**
	@brief * operator for getting the value of the iterator
	@return FASTAEntry with specific candidate
	@throw Exception::InvalidIterator if iterator has not been initialized
	*/
	virtual FASTAEntry operator*() throw (Exception::InvalidIterator);
	
	/**
	@brief opperator ++ for postincrement
	@return Reference to PepIterator
	@throw Exception::InvalidIterator if iterator has not been initialized
	*/
	virtual PepIterator & operator++() throw (Exception::InvalidIterator);
	
	/**
	@brief opperator ++ for preincrement
	@return pointer to PepIterator
	@throw Exception::InvalidIterator if iterator has not been initialized
	*/
	virtual PepIterator * operator++(int i) throw (Exception::InvalidIterator);

	/**
	@brief setter for fasta file
	@param f String with fasta file location
	@throw Exception::FileNotFound if file could not be found
	*/
	virtual void setFastaFile (const String & f) throw (Exception::FileNotFound);
	
	/**
	@brief getter for FASTA file
	@return String with file location
	*/
	virtual String getFastaFile();

	/**
	@brief setter for tolerance
	@throw Exception::NotImplemented because its not available for tryptic iterator
	*/
	virtual void setTolerance (float) throw (Exception::InvalidValue, Exception::NotImplemented)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	};
	
	/**
	@brief getter for tolerance
	@return tolerance
	@throw Exception::NotImplemented because its not available for tryptic iterator
	*/
	virtual float getTolerance () throw (Exception::InvalidValue, Exception::NotImplemented)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	};

	/**
	@brief setter for spectrum
	@throw Exception::NotImplemented because its not available for tryptic iterator
	*/
	virtual void setSpectrum (const std::vector<float> & ) throw (Exception::InvalidValue, Exception::NotImplemented)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	};

	/**
	@brief getter for spectrum
	@return the used spectrum
	@throw Exception::NotImplemented because its not available for tryptic iterator
	*/
	virtual const std::vector<float> & getSpectrum () throw (Exception::InvalidValue, Exception::NotImplemented)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	};

	/**
	@brief initializing iterator
	@return true if everything was ok
	*/
	virtual bool begin () throw (Exception::InvalidIterator);

	/**
	@brief indicates whether iterator is at end
	@return true if iterator is at end
	@see hasNext
	*/
	virtual bool isAtEnd ();

	/**
	@brief indicated if a digesting enzyme will cut at this position
	@return true if digenting enzym cuts the sequence
	*/
	virtual bool isDigestingEnd (char aa1, char aa2);
	
	/**
	@brief needed by FactoryProduct
	@return const string name of class
	*/
	static const std::string getName()
	{
		return "TrypticIterator";
	}

	/**
	@brief needed by FactoryProduct
	@return poiter to new object
	*/
	static PepIterator* create()
	{
		return new TrypticIterator;
	}
	
	
	
	protected:
		/**
		@brief getting the next candidate
		@return string with next sequence
		*/
		virtual std::string next_ ();
	
		/**
		@brief indicates if there will be a next element
		@return true if iterator has more elements
		*/
		bool hasNext_();
	
		/**
		@brief finds the next starting position where a digesting enzyme will cut the sequence
		*/
		void goToNextAA_ ();

		String f_file_; ///< fasta file location

		std::string actual_pep_; ///< actual peptide

		bool is_at_end_; ///< indicates if iterator is at end

		PepIterator* f_iterator_; ///< FastaIterator

		FASTAEntry f_entry_; ///< actual fasta entry

		unsigned int b_,e_; ///< to ints representing a position within the actual string (b = begin, e = end)

	
	
};

}
#endif //OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATOR_H
