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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/BigString.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>

namespace OpenMS {

typedef std::pair<String,String> FASTAEntry;

BigString::BigString () :
	big_string_("$"),
	separator_('$'),
	count_(1),
	len_(1)
{
	sep_indices_.push_back(0);
	FASTA_header_.push_back("");
}

BigString::BigString (const BigString & bs):
	big_string_(bs.big_string_),
	separator_(bs.separator_),
	count_(bs.count_),
	len_(bs.len_),
	sep_indices_(bs.sep_indices_),
	FASTA_header_(bs.FASTA_header_)
{

}

BigString::~BigString()
{

}

void BigString::add (FASTAEntry const & new_entry)
{
	big_string_+=new_entry.second;
	big_string_+=separator_;
	++count_;
	len_ += new_entry.second.length() +1;
	sep_indices_.push_back(len_-1);
	FASTA_header_.push_back(new_entry.first);
}

void BigString::setSeparator (const char sep)
{
	separator_ = sep;
}

char BigString::getSeparator ()
{
	return (separator_);
}

Size BigString::size ()
{
	return (count_);
}

Size BigString::length ()
{
	return (len_);
}

void BigString::getPeptide(FASTAEntry& entry, Size start, Size length)
{
	Size index_start = getIndex_(start);
	if (index_start != getIndex_(start + length))
	{
		throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "desired peptide is part of 2 fasta entries","");
	}
	entry.first = FASTA_header_[index_start];
	entry.second = big_string_.substr(start, length);
	return;
}

const String& BigString::getBigString() const
{
	return big_string_;
}

Size BigString::getIndex_(Size index, Size start, Size end)
{
	if (end - start <= 1)
	{
		if (sep_indices_[start] >= index)
		{
			return start;
		}
		else
		{
			return start + 1;
		}
	}
	Size half =(Size) ((end-start)/2)+start;

	if (index > sep_indices_[half])
	{
		return getIndex_(index, half, end);
	}
	else if (index < sep_indices_[half])
	{
		return getIndex_(index, start, half);
	}
	else
	{
		return half;
	}
}

Size BigString::getIndex_(Size index)
{
	return getIndex_(index, 0, sep_indices_.size());
}

}
