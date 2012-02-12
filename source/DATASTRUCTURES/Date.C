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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Date.h>

#include <QtCore/QDate>

#include <ctime>

using namespace std;

namespace OpenMS
{
	Date::Date():
		QDate()
	{
		
	}

	Date::Date(const Date& date): 
		QDate(date)
	{
	}

	Date::Date(const QDate& date): 
		QDate(date)
	{
	}

	Date& Date::operator= (const Date& source)
	{
	  if (&source == this)
	  { 
	  	return *this;
	  }
	  QDate::operator=(source);
	  
	  return *this;		
	}

	void Date::set(const String& date)
	{
		clear();
		
		//check for format (german/english)
		if (date.has('.'))
		{
			QDate::operator=(QDate::fromString(date.c_str(), "dd.MM.yyyy"));
		}
		else if (date.has('/'))
		{
			QDate::operator=(QDate::fromString(date.c_str(), "MM/dd/yyyy"));
		}
		else if (date.has('-'))
		{
			QDate::operator=(QDate::fromString(date.c_str(), "yyyy-MM-dd"));
		}

		if (!isValid())
		{
			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__, date, "Is no valid german, english or iso date");
		}
	}
	
	void Date::set(UInt month, UInt day, UInt year)
	{
		if (!setDate(year, month, day))
		{	
			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__, String(year) + "-" + String(month) + "-" + String(day), "Invalid date");
		}
	}

	Date Date::today()
	{
		return QDate::currentDate();
	}

	String Date::get() const
	{
		if (QDate::isValid())
		{
			return toString("yyyy-MM-dd");
		}
		return "0000-00-00";
	}
	
	void Date::get(UInt& month, UInt& day, UInt& year) const
	{
		day = QDate::day();
		month = QDate::month();
		year = QDate::year();
	}

	void Date::clear()
	{
		QDate::operator=(QDate());
	}
	
} // namespace OpenMS


