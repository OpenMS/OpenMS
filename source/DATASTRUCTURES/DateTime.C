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
// $Maintainer: Sandro Andreotti $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <iostream>
#include <ctime>

using namespace std;

namespace OpenMS
{
	DateTime::DateTime(): 
		QDateTime()
	{
		
	}

	DateTime::DateTime(const DateTime& date): 
		QDateTime(date)
	{
		
	}

	DateTime::DateTime(const QDateTime& date): 
		QDateTime(date)
	{
		
	}

	DateTime& DateTime::operator= (const DateTime& source)
	{
	  if (&source == this)
	  { 
	  	return *this;
	  }
	  
	  QDateTime::operator=(source);
	  
	  return *this;		
	}

	void DateTime::set(const String& date)
	{
		clear();
		
		if (date.has('.'))
		{
			QDateTime::operator=(QDateTime::fromString(date.c_str(), "dd.MM.yyyy hh:mm:ss"));
		}
		else if (date.has('/'))
		{
			QDateTime::operator=(QDateTime::fromString(date.c_str(), "MM/dd/yyyy hh:mm:ss"));
		}
		else if (date.has('-'))
		{
			if (date.has('T'))
			{
				QDateTime::operator=(QDateTime::fromString(date.c_str(), "yyyy-MM-ddThh:mm:ss"));
			}
			else if (date.has('Z'))
			{
				QDateTime::operator=(QDateTime::fromString(date.c_str(), "yyyy-MM-ddZ"));
			}
			else if (date.has('+'))
			{
				QDateTime::operator=(QDateTime::fromString(date.c_str(), "yyyy-MM-dd+hh:mm"));
			}
			else
			{
				QDateTime::operator=(QDateTime::fromString(date.c_str(), "yyyy-MM-dd hh:mm:ss"));
			}
		}
		
		if (!QDateTime::isValid())
		{
    	throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Invalid date time string");			
		}
	}
	
	void DateTime::set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second)
	{
		QDateTime::setDate(QDate(year, month, day));
		QDateTime::setTime(QTime(hour, minute, second));

		if (!QDateTime::isValid())
		{
			String date_time = String(year) + "-" + String(month) + "-" + String(day) 
				+ " " + String(hour) + ":" + String(minute) + ":" + String(second);
			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,date_time,"Invalid date time");
		}
	}

	DateTime DateTime::now()
	{
		return QDateTime::currentDateTime();
	}

	String DateTime::get() const
	{
		if (QDateTime::isValid())
		{
			return QDateTime::toString("yyyy-MM-dd hh:mm:ss");
		}
		return "0000-00-00 00:00:00";
	}
	
	void DateTime::get(UInt& month, UInt& day, UInt& year,
										UInt& hour, UInt& minute, UInt& second) const
	{
		const QDate& temp_date = QDateTime::date();
		const QTime& temp_time = QDateTime::time();
		
		year = temp_date.year();
		month = temp_date.month();
		day = temp_date.day();
		hour = temp_time.hour();
		minute = temp_time.minute();
		second = temp_time.second();		
	}

	void DateTime::clear()
	{
		QDateTime::operator=(QDateTime());
	}

	void DateTime::setDate(const String& date)
	{
		QDate temp_date;
		
		if (date.has('-'))
		{
			temp_date = QDate::fromString(date.c_str(), "yyyy-MM-dd");
		}
		else if (date.has('.'))
		{			
			temp_date = QDate::fromString(date.c_str(), "dd-MM-yyyy");
		}
		else if (date.has('/'))
		{			
			temp_date = QDate::fromString(date.c_str(), "MM/dd/yyyy");
		}
		else
		{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Could not set date");
		}
		if (!temp_date.isValid())
		{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Could not set date");
		}
		
		QDateTime::setDate(temp_date);
	}
				
	void DateTime::setTime(const String& time)
	{
		QTime temp_time;
		
		temp_time = QTime::fromString(time.c_str(), "hh:mm:ss");
		if (!temp_time.isValid())
		{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, time, "Could not set time");
		}
		
		QDateTime::setTime(temp_time);				
	}
				
	void DateTime::setDate(UInt month, UInt day, UInt year)
	{
		QDate temp_date;
		
		if (!temp_date.setDate(year, month, day))
		{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String(year) + "-" + String(month) + "-" + String(day), "Could not set date");
		}
		
		QDateTime::setDate(temp_date);		
	}
		
	void DateTime::setTime(UInt hour, UInt minute, UInt second)
	{
		QTime temp_time;
		
		if (!temp_time.setHMS(hour, minute, second))
		{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String(hour) + ":" + String(minute) + ":" + String(second), "Could not set time");
		}
		QDateTime::setTime(temp_time);
	}
	
	void DateTime::getDate(UInt& month, UInt& day, UInt& year) const
	{
		const QDate& temp_date = QDateTime::date();
		
		month = temp_date.month();
		day = temp_date.day();
		year = temp_date.year();
	}

	String DateTime::getDate() const
	{
		if (QDateTime::isValid())
		{
			return QDateTime::date().toString("yyyy-MM-dd");
		}
		return "0000-00-00";
	}

	void DateTime::getTime(UInt& hour, UInt& minute, UInt& second) const
	{
		const QTime& temp_time = QDateTime::time();
		
		hour = temp_time.hour();
		minute = temp_time.minute();
		second = temp_time.second();
	}

	String DateTime::getTime() const
	{
		if (QDateTime::isValid())
		{
			return QDateTime::time().toString("hh:mm:ss");
		}
		return "00:00:00";
	}
			
} // namespace OpenMS


