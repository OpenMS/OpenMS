// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

FastaIterator::FastaIterator() :  PepIterator()
{
	actual_seq_ = "";
	is_at_end_ = false;
	const String ff = "";
	fasta_file_= ff;
	input_file_ = new std::ifstream();
}

FastaIterator::~FastaIterator()
{
	
}

FastaIterator::FastaIterator(const FastaIterator & source) : PepIterator(source)
{
	is_at_end_ = (source.is_at_end_);
	input_file_ = (source.input_file_);
	fasta_file_ = (source.fasta_file_);
	actual_seq_ =(source.actual_seq_);
	header_ = (source.header_);
	last_header_ = (source.last_header_);
}

PepIterator * FastaIterator::operator++(int)
{
	if (last_header_=="")
	{
		throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	PepIterator * old = new FastaIterator (*this);
	actual_seq_ = next_();
	return old;
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
	fs.open(f.c_str());
	if (!fs.is_open())
	{
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, f);
	}
	fasta_file_ = f;
}

String FastaIterator::getFastaFile()
{
	return (fasta_file_);
}

std::string FastaIterator::next_()
{
	if (input_file_->eof())
	{
		is_at_end_ = true;
		return ("");
	}
	std::string line;
	std::getline(*input_file_, line);
	if (line[0] == '>' || input_file_->eof())
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
	input_file_->open(fasta_file_.c_str());
	
	if (*input_file_)
	{
		std::string line;
		std::getline(*input_file_,line);
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
