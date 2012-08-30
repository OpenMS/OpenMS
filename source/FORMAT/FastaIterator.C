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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/FastaIterator.h>
#include <iostream>


namespace OpenMS
{

typedef std::pair <String, String> FASTAEntry;

FastaIterator::FastaIterator() :  
  PepIterator(),
  actual_seq_(),
  is_at_end_(false),
  fasta_file_(),
  input_file_()
{
}

FastaIterator::~FastaIterator()
{
	
}

// not implemented (since copying this stuff will invalidate the istream
//FastaIterator::FastaIterator(const FastaIterator & source) : PepIterator(source)
//{
//}

// not implemented (since copying this stuff will invalidate the istream
//FastaIterator& operator=(const FastaIterator &);


PepIterator * FastaIterator::operator++(int)
{ // this operator requires copying, which we cannot support!
	throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
}

FASTAEntry FastaIterator::operator*() 
{
	if (last_header_=="")
	{
		throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	return FASTAEntry (last_header_,actual_seq_);
}

PepIterator & FastaIterator::operator++()
{
	if (last_header_=="")
	{
		throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	actual_seq_ = next_();
	return *this;	
}

void FastaIterator::setFastaFile (const String & f)
{
	std::fstream fs;
	fs.open(f.c_str(), std::fstream::in);
	if (!fs.is_open())
	{
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, f);
	}
  fs.close();
	fasta_file_ = f;
}

String FastaIterator::getFastaFile()
{
	return (fasta_file_);
}

std::string FastaIterator::next_()
{
	if (input_file_.eof())
	{
		is_at_end_ = true;
    input_file_.close();
		return ("");
	}
  is_at_end_ = false;
	std::string line;
	std::getline(input_file_, line);
	if (line[0] == '>' || input_file_.eof())
	{
		last_header_ = header_;
		header_ = line;
		return ("");
	}
	return (std::string(line)+next_());
}
	
bool FastaIterator::begin()
{
	if (fasta_file_=="")
	{
		throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	input_file_.open(fasta_file_.c_str(), std::fstream::in);
	
	if (input_file_)
	{
		std::string line;
		std::getline(input_file_,line);
		header_ = line;
		last_header_ = line;
		actual_seq_ = next_();
		return (true);
	}
	
	return (false);
	
}
	
bool FastaIterator::isAtEnd ()
{
	return (is_at_end_);	
}

} //namespace OpenMS
