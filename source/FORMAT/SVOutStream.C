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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <limits>

using namespace std;

namespace OpenMS
{

SVOutStream::SVOutStream(ostream& out, const String& sep,
												 const String& replacement,
												 String::QuotingMethod quoting):
	ostream(out.rdbuf()), sep_(sep), replacement_(replacement),	nan_("nan"), 
	inf_("inf"), quoting_(quoting), modify_strings_(true), newline_(true)
{
	// use high decimal precision (appropriate for double):
	precision(std::numeric_limits<double>::digits10);
}


SVOutStream& SVOutStream::operator<<(String str)
{
	if (str.find('\n') != String::npos)
	{
		throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "argument must not contain newline characters");
	}

	if (!newline_) 
	{
		(ostream&)*this << sep_;
	}
	else
	{
		newline_ = false;
	}
	
	if (!modify_strings_)
	{
		(ostream&)*this << str;
	}
	else if (quoting_ != String::NONE)
	{
		(ostream&)*this << str.quote('"', quoting_);
	}
	else 
	{
		(ostream&)*this << str.substitute(sep_, replacement_);
	}
	return *this;
}


SVOutStream& SVOutStream::operator<<(const char* c_str)
{
	return operator<<(String(c_str));
}


SVOutStream& SVOutStream::operator<<(const char c)
{
	return operator<<(String(c));
}


SVOutStream& SVOutStream::operator<<(const string& str)
{
	return operator<<((String&)str);
}


SVOutStream& SVOutStream::operator<<(ostream& (*fp)(ostream&))
{
	ostream& (* const endlPointer) (ostream&) = &endl;
	if (fp == endlPointer)	newline_ = true;
	(ostream&)*this << fp;
	return *this;
}


SVOutStream& SVOutStream::write(const String& str)
{
	ostream::write(str.c_str(), str.size());
	return *this;
}


bool SVOutStream::modifyStrings(bool modify)
{
	bool old = modify_strings_;
	modify_strings_ = modify;
	return old;
}

} // namespace OpenMS


