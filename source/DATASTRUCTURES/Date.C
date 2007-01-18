// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Date.h>

#include <time.h>

using namespace std;

namespace OpenMS
{
	Date::Date(): 
		day_(0),
		month_(0),
		year_(0)
	{
		
	}

	Date::Date(const Date& date): 
		day_(date.day_),
		month_(date.month_),
		year_(date.year_)
	{
		
	}

	Date::~Date()
	{
		
	}
	
	Date& Date::operator= (const Date& source)
	{
	  if (&source == this) return *this;
	  
		day_ = source.day_;
		month_ = source.month_;
		year_ = source.year_;
	  
	  return *this;		
	}

	bool Date::operator == (const Date& rhs) const	
	{
		return year_==rhs.year_ && month_==rhs.month_ && day_==rhs.day_;
	}

	bool Date::operator != (const Date& rhs) const	
	{
		return !(operator==(rhs));
	}
	
	void Date::set(const String& date) throw (Exception::ParseError)
	{
		std::vector<String> split;
		stringstream ss;
		String d, m, y;
		UnsignedInt day, month, year;
		
		//check for format (german/english)
		if (date.has('.'))
		{
			date.split('.',split);

			//check for right number of parts
			if (split.size()!=3)
			{
				throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,date, "Is no valid german date");
			}
			d=split[0];
			m=split[1];
			y=split[2];
			
		}
		else if (date.has('/'))
		{
			date.split('/',split);

			//check for right number of parts
			if (split.size()!=3)
			{
				throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,date, "Is no valid german date");
			}
			d=split[1];
			m=split[0];
			y=split[2];
		}
		else if (date.has('-'))
		{
			date.split('-',split);

			//check for right number of parts
			if (split.size()!=3)
			{
				throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,date, "Is no valid iso date");
			}
			d=split[2];
			m=split[1];
			y=split[0];
		}
		else
		{
			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__, date, "Is no valid german, english or iso date");
		}
		
		//day
    ss << d;
    if (!(ss >> day))
    {
    	throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Could not convert day to a number");
    } 
    
    //month
    ss.clear();
    ss << m;
    if (!(ss >> month))
    {
    	throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Could not convert month to a number");
    }
    
    //year
    ss.clear();
    ss << y;
    if (!(ss >> year))
    {
    	throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Could not convert year to a number");
    }
		
		set(month,day,year);
	}
	
	void Date::set(UnsignedInt month, UnsignedInt day, UnsignedInt year) throw (Exception::ParseError)
	{
		//set month lengths
		UnsignedInt month_length[13] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    if ( isLeapYear(year) )
    {
    	month_length[2] = 29;
    }

		//month correct
		if (month==0 || month>12)
		{
			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,String(month),"Invalid month");
		}
		
		//day correct
		if (day==0 || day>month_length[month])
		{
			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__, String(day), "Invalid day");
		}

		day_ = day;
		month_ = month;
		year_ = year;				
	}

	void Date::today()
	{
		time_t sec = time(NULL);
		struct tm* tmp =  localtime(&sec);
		set(tmp->tm_mon+1, tmp->tm_mday, 1900+tmp->tm_year);		
	}

	string Date::now()
	{
		char out[60]; 
		//get time in sec
		time_t sec = time(NULL);
		struct tm* ptr = localtime(&sec);
		//struct tm* tmp =  ;
		strftime( out, 100, "%Y-%m-%d %H:%M:%S" , ptr );
		
		return out;
	}
	
	void Date::get(String& date) const
	{
		date = String(year_).fillLeft('0',4)+"-"+String(month_).fillLeft('0',2)+"-"+String(day_).fillLeft('0',2);
	}
	
	void Date::get(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year) const
	{
		day = day_;
		month = month_;
		year = year_;
	}

	void Date::clear()
	{
		day_ = 0;
		month_ = 0;
		year_ = 0;
	}
	
	bool Date::isLeapYear(UnsignedInt year) const
	{
 		if ( (year%4) != 0 ) return false;
		if ( (year%400) == 0 ) return true;
		if ( (year%100) == 0 ) return false;
		return true;
	}

} // namespace OpenMS


