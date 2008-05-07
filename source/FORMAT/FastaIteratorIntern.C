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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/FastaIteratorIntern.h>
#include <OpenMS/FORMAT/FASTAFile.h>

namespace OpenMS {

	typedef std::pair<String,String> FASTAEntry;
	
	FastaIteratorIntern::FastaIteratorIntern() :
		fasta_file_("")
	{
		
	}
	
	FastaIteratorIntern::~FastaIteratorIntern()
	{
		
	}

	FastaIteratorIntern::FastaIteratorIntern(const FastaIteratorIntern & source) :
		PepIterator (source),
		fasta_file_(source.fasta_file_),
		entrys_(source.entrys_),
		it_(source.it_)
	{
		
	}

	FASTAEntry FastaIteratorIntern::operator*() throw (Exception::InvalidIterator)
	{
		if (fasta_file_=="")
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		return *it_;
	}
	
	PepIterator & FastaIteratorIntern::operator++() throw (Exception::InvalidIterator)
	{
		if (fasta_file_=="")
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		++it_;
		return *this;
	}
	
	PepIterator * FastaIteratorIntern::operator++(int) throw (Exception::InvalidIterator)
	{
		if (fasta_file_=="")
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		PepIterator * old = new FastaIteratorIntern (*this);
		++it_;
		return old;
	}
	
	void FastaIteratorIntern::setFastaFile (const String & f) throw (Exception::FileNotFound, Exception::ParseError)
	{
		
		FASTAFile ffile;
		ffile.load (f,entrys_);
		fasta_file_ = f;
		it_ = entrys_.begin();
	}

	String FastaIteratorIntern::getFastaFile ()
	{
		return (fasta_file_);
	}
	
	bool FastaIteratorIntern::begin () throw (Exception::InvalidIterator)
	{
		if (fasta_file_=="")
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		return true;
	}
	
	bool FastaIteratorIntern::isAtEnd ()
	{
		return (it_ == entrys_.end());
	}

}
