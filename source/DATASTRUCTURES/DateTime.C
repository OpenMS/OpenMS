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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <iostream>
#include <time.h>

using namespace std;

namespace OpenMS
{
	DateTime::DateTime(): 
		Date(),
		hour_(0),
		minute_(0),
		second_(0)
	{
		
	}

	DateTime::DateTime(const DateTime& date): 
		Date(date),
		hour_(date.hour_),
		minute_(date.minute_),
		second_(date.second_)
	{
		
	}

	DateTime::~DateTime()
	{
		
	}
	
	DateTime& DateTime::operator= (const DateTime& source)
	{
	  if (&source == this) return *this;
	  
		day_ = source.day_;
		month_ = source.month_;
		year_ = source.year_;
		hour_ = source.hour_;
		minute_ = source.minute_;
		second_ = source.second_;
	  
	  return *this;		
	}

	bool DateTime::operator == (const DateTime& rhs) const	
	{
		return year_==rhs.year_ && month_==rhs.month_ && day_==rhs.day_
					&& hour_ == rhs.hour_ && minute_ == rhs.minute_ 
					&& second_ == rhs.second_;
	}

	bool DateTime::operator != (const DateTime& rhs) const	
	{
		return !(operator==(rhs));
	}
	
	void DateTime::set(const String& date) throw (Exception::ParseError)
	{
		std::vector<String> split;
		
		//separate date and time
		date.split(' ',split);
		
		// if no time is given only set the date
		if (split.size() == 0)
		{
			Date::set(date);
		}
		else if (split.size() == 2)
		{
			// set the date
			Date::set(split[0].trim());
			
			// set the time
    	setTime(split[1].trim()); 
		}
		else
		{
    		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Invalid date time string");			
		}
	}
	
	void DateTime::set(UnsignedInt month, UnsignedInt day, UnsignedInt year,
							 			 UnsignedInt hour, UnsignedInt minute, UnsignedInt second) throw (Exception::ParseError)
	{
		Date::set(month, day, year);

		if (hour > 23)
		{
			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,String(hour),"Invalid hour");
		}
		if (minute > 59)
		{
			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,String(hour),"Invalid minute");
		}
		// second can also be 61 to accout for leap seconds (like in tm_struct)
		if (second > 61)
		{
			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,String(hour),"Invalid second");
		}
		setTime(hour, minute, second);
	}

	void DateTime::now()
	{
		time_t sec = time(NULL);
		struct tm* tmp =  localtime(&sec);
		setDate(tmp->tm_mon+1, tmp->tm_mday, 1900+tmp->tm_year);
		setTime(tmp->tm_hour, tmp->tm_min, tmp->tm_sec);
	}

	void DateTime::get(String& date) const
	{
		String partial_date;
		Date::get(partial_date);
		
		date = partial_date + " ";
		if (hour_ < 10)
		{
			date = date + "0";
		}
		date = date + String(hour_) + ":";
		if (minute_ < 10)
		{
			date = date + "0";
		}
		date = date + String(minute_) + ":";
		if (second_ < 10)
		{
			date = date + "0";			
		}
		date = date + String(second_);
	}
	
	void DateTime::get(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year,
										UnsignedInt& hour, UnsignedInt& minute, UnsignedInt& second) const
	{
		Date::get(month, day, year);
		hour = hour_;
		minute = minute_;
		second = second_;		
	}

	void DateTime::clear()
	{
		Date::clear();
			
		hour_ = 0;
		minute_ = 0;
		second_ = 0;
	}

	void DateTime::setDate(const String& date) throw (Exception::ParseError)
	{
		Date::set(date);
	}
				
	void DateTime::setTime(const String& date) throw (Exception::ParseError)
	{
		vector<String> split;
		String h, m, s;
		UnsignedInt hour, minute, second;
		stringstream ss;	
		
		date.split(':', split);
		if (split.size() != 3)
		{
			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__, date, "The time information is wrong");			
		}
		h = split[0];
		m = split[1];
		s = split[2];
		
		//hour
    ss << h;
	  if (!(ss >> hour))
  	{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Could not convert hour to a number");
  	} 

    ss.clear();
		//minute
    ss << m;
	  if (!(ss >> minute))
  	{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Could not convert minute to a number");
  	} 

    ss.clear();
		//second
    ss << s;
	  if (!(ss >> second))
  	{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Could not convert second to a number");
  	}
   	setTime(hour, minute, second); 		
	}
				
	void DateTime::setDate(UnsignedInt month, UnsignedInt day, UnsignedInt year) throw (Exception::ParseError)
	{
		Date::set(month, day, year);
	}
		
	void DateTime::setTime(UnsignedInt hour, UnsignedInt minute, UnsignedInt second) throw (Exception::ParseError)
	{
		hour_ = hour;
		minute_ = minute;
		second_ = second;
	}
	
	void DateTime::getDate(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year) const
	{
		month = month_;
		day = day_;
		year = year_;
	}

	void DateTime::getDate(String& date) const
	{
		Date::get(date);
	}

	void DateTime::getTime(UnsignedInt& hour, UnsignedInt& minute, UnsignedInt& second) const
	{
		hour = hour_;
		minute = minute_;
		second = second_;
	}

	void DateTime::getTime(String& time) const
	{
		time.clear();
		time = time + String(hour_) + ":" + String(minute_) + ":" + String(second_);		
	}
			
} // namespace OpenMS


