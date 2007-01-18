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

#include <OpenMS/DATASTRUCTURES/String.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <limits>
#include <algorithm>

using namespace std;

namespace OpenMS
{
	const String String::EMPTY;
	
	const SignedInt String::NPOS(std::numeric_limits<SignedInt>::min());

	String::String()
		:	string()
	{
	}

	String::String(const string& s)
		:	string(s)
	{
	}

	String::String(const char* s)
		:	string(s)
	{
	}

	String::String(const char* s, SizeType length)
	{
		SizeType count = 0;
		while(count < length && *(s+count)!=0)
		{
			*this += *(s+count);
			++count;
		}
	}

	String::String(const char c)
		:	string(1,c)
	{
	}

	String::String(size_t len, char c)
		:	string(len, c)
	{
	}

	String::String(int i)
		: string()
	{
		stringstream s;
		s << i;
		string::operator=(s.str());
	}

	String::String(unsigned int i)
		: string()
	{
		stringstream s;
		s << i;
		string::operator=(s.str());
	}

	String::String(short int i)
		: string()
	{
		stringstream s;
		s << i;
		string::operator=(s.str());
	}

	String::String(short unsigned int i)
		: string()
	{
		stringstream s;
		s << i;
		string::operator=(s.str());
	}

	String::String(long int i)
		: string()
	{
		stringstream s;
		s << i;
		string::operator=(s.str());
	}

	String::String(long unsigned int i)
		: string()
	{
		stringstream s;
		s << i;
		string::operator=(s.str());
	}

	String::String(long long unsigned int i)
		: string()
	{
		stringstream s;
		s << i;
		string::operator=(s.str());
	}

	String::String(float f)
		: string()
	{
		stringstream s;
		s.precision(7);
		s << f;
		string::operator=(s.str());
	}

	String::String(double d)
		: string()
	{
		stringstream s;
		s.precision(10);
		s << d;
		string::operator=(s.str());
	}

	String::String(long double d)
		: string()
	{
		stringstream s;
		s.precision(16);
		s << d;
		string::operator=(s.str());
	}
	
	String::String(double d, UnsignedInt size)
		: string()
	{
		stringstream s;
		//reserve one space for the minus sign
		SignedInt sign=0;
		if (d<0) sign=1;
		d = fabs(d);
		
		if (d<pow(10.0,SignedInt(size-sign-2)))
		{
			s.precision(10);
			if (sign==1) s << "-";
			s << d;
		}
		else
		{
			UnsignedInt exp=0;
			while(d>pow(10.0,SignedInt(size-sign-4)))
			{
				d/=10;
				++exp;
			}
			d = SignedInt(d)/10.0;
			exp+=1;
			if (sign==1) s << "-";
			s << d<<"e";
			if (exp<10) s << "0";
			s <<exp; 
		}
		string::operator=(s.str().substr(0,size));		
	}

	String& String::fillLeft(char c, UnsignedInt size)
	{
		if (string::size()<size)
		{
			string::operator=(string(size-string::size(),c)+*this);
		}
		return *this;
	}

	String& String::fillRight(char c, UnsignedInt size)
	{
		if (string::size()<size)
		{
			string::operator=(*this+string(size-string::size(),c));
		}		
		return *this;
	}
	
	String::String(const DataValue& d)
		: string()
	{
		string::operator=(d.toString());
	}

	bool String::hasPrefix(const String& string) const
	{
		if (string.size() > size())
		{
			return false;
		}
		if (string.size() == 0)
		{
			return true;
		}
		return (compare(0, string.size(), string) == 0);
	}

	bool String::hasSuffix(const String& string) const
	{
		if (string.size() > size())
		{
			return false;
		}
		if (string.size() == 0)
		{
			return true;
		}
		return (compare(size()-string.size(), string.size(), string) == 0);
	}

	bool String::hasSubstring(const String& string) const
	{
		if (string.size() > size())
		{
			return false;
		}
		if (string.size() == 0)
		{
			return true;
		}
		for (Position i=0;i!=size();++i)
		{
			if (compare(i, string.size(), string) == 0)
			{
				return true;
			}
		}
		return false;
	}

	bool String::has(Byte byte) const
	{
		for (Position i=0;i!=size();++i)
		{
			if ((*this)[i] == byte)
			{
				return true;
			}
		}
		return false;
	}

	String String::prefix(SizeType length) const
		throw(Exception::IndexOverflow)
	{
		if (length > size())
		{
			throw(Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, length, size()));
		}
		return substr(0, length);
	}

	String String::suffix(SizeType length) const
		throw(Exception::IndexOverflow)
	{
		if (length > size())
		{
			throw(Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, length, size()));
		}
		return substr(size()-length, length);
	}

	String String::prefix(SignedInt length) const
		throw(Exception::IndexUnderflow, Exception::IndexOverflow)
	{
		if (length < 0)
		{
			throw(Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, length,0));
		}
		if (length > SignedInt(size()))
		{
			throw(Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, length, size()));
		}
		return substr(0, length);
	}

	String String::suffix(SignedInt length) const
		throw(Exception::IndexUnderflow, Exception::IndexOverflow)
	{
		if (length < 0)
		{
			throw(Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, length,0));
		}
		if (length > SignedInt(size()))
		{
			throw(Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, length, size()));
		}
		return substr(size()-length, length);
	}

	String String::prefix(char delim) const
		throw(Exception::ElementNotFound<char>)
	{
		ConstIterator it=begin();
		while (it!=end() && *it!=delim)
		{
			++it;
		}
		//char not found
		if (it==end())
		{
			throw(Exception::ElementNotFound<char>(__FILE__, __LINE__, __PRETTY_FUNCTION__, delim));
		}
		return String(begin(), it);
	}

	String String::suffix(char delim) const
		throw(Exception::ElementNotFound<char>)
	{
		ConstIterator it=end();
		--it;
		while (it!=--(begin()) && *it!=delim)
		{
			--it;
		}
		//char not found
		if (it==--(begin()))
		{
			throw(Exception::ElementNotFound<char>(__FILE__, __LINE__, __PRETTY_FUNCTION__, delim));
		}
		++it;
		return String(it, end());
	}

	String String::substr(SignedInt start, SignedInt n) const
	{
		if (start>=0 && (n>=0 || n==NPOS))
		{
			return string::substr(start,n);
		}
		
		ConstIterator begin;
		if (start<0)
		{
			begin = this->end()+start;
		}
		else
		{
			begin = this->begin()+start;
		}
		
		ConstIterator end;
		if (n==NPOS)
		{
			end = this->end();
		}
		else if (n<0)
		{
			end = this->end()+n;
			if (end<begin)
			{
				end = begin;
			}
		}
		else
		{
			end = begin + n;
		}		
		
		return String(begin,end);
		
	}

	String& String::trim()
	{
		//search for the begin of truncated string
		iterator begin = this->begin();
		while (begin!=end() && (*begin==' ' || *begin=='\t' || *begin=='\n'  || *begin=='\r' ))
		{
			++begin;
		}
		
		//all characters are whitespaces
		if (begin==end())
		{
			string::clear();
			return *this;
		}
		
		//search for the end of truncated string
		Iterator end=this->end();
		end--;
		while (end!=begin && (*end==' ' || *end=='\n' || *end=='\t' || *end=='\r' ))
		{
			--end;
		}
		++end;

		//no characters are whitespaces
		if (begin==this->begin() && end==this->end())
		{
			return *this;
		}
		
		string::operator=(string(begin,end));
		return *this;
	}

	String String::random(Size length)
	{
		srand(time(0));
		String tmp(length, '.');
		Size random;
		for (Position i = 0 ; i < length; ++i)
		{
			random = (SizeType)floor(((double)rand()/(double(RAND_MAX)+1)) * 62.0);
			if (random < 10)
			{
				tmp[i] = (char)(random +48);
			}
			else if (random < 36)
			{
				tmp[i] = (char)(random +55);
			}
			else
			{
				tmp[i] = (char)(random +61);
			}
		}
		return tmp;
	}

	String& String::reverse()
	{
		String tmp = *this;
		for (Position i=0;i!=size();++i)
		{
			(*this)[i] = tmp[size()-1-i];
		}
		return *this;
	}

	bool String::split(char splitter, std::vector<String>& substrings) const
	{
		SignedInt parts = count(this->begin(),this->end(),splitter);
		
		// no splitter found
		if (parts == 0)
		{
			substrings.clear();
			return false;
		}
		
		// splitter(s) found
		substrings.resize(parts+1);		
		
		ConstIterator end = this->begin();
		ConstIterator begin = this->begin();		
		
		parts=0;
		
		for (; end != this->end(); ++end)
		{
			if (*end == splitter)
			{
				substrings[parts++] = String(begin,end);
				begin = end+1;
			}
		}
		substrings[parts++] = String(begin,end);
		
		return true;
	}

	void String::implode(vector<String>::const_iterator first, vector<String>::const_iterator last, const string& glue)
	{
		//empty container
		if (first==last)
		{
			string::clear();
			return;
		}
		
		string::operator=(*first);
		for (vector<String>::const_iterator it = ++first; it != last; ++it)
		{
			string::operator+=( glue + (*it));
		}
	}


	int String::toInt() const throw(Exception::ConversionError)
	{
    std::stringstream ss(c_str());
    int ret;
    if (!(ss >> ret)) throw(Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Could not convert string ")+*this+" to an int"));
    return ret;
    // long int strtol(const char *nptr, char **endptr, int base);    
    // return atoi(c_str());
	}

	// long long int strtoll(const char *nptr, char **endptr, int base); // probably not easily portable


	float String::toFloat() const throw(Exception::ConversionError)
	{
    std::stringstream ss(c_str());
    float ret;
    if (!(ss >> ret)) throw(Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,  String("Could not convert string ")+*this+" to a float"));
    return ret;    
		// float strtof(const char *nptr, char **endptr);
		// return (float)atof(c_str());
	}

	double String::toDouble() const throw(Exception::ConversionError)
	{
    std::stringstream ss(c_str());
    double ret;
    if (!(ss >> ret)) throw(Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,  String("Could not convert string ")+*this+" to a double"));
    return ret;    
		// double strtod(const char *nptr, char **endptr);
		//return atof(c_str());
	}

	String& String::toUpper()
	{
		std::transform(this->begin(), this->end(), this->begin(), (int(*)(int)) toupper);
		return *this;
	}

	String& String::firstToUpper()
	{
		if (this->size()!=0)
		{
			(*this)[0] = toupper ((*this)[0]);
		}
		return *this;
	}

	String& String::toLower()
	{
		std::transform(this->begin(), this->end(), this->begin(), (int(*)(int)) tolower);
		return *this;
	}
	
	String& String::substitute(char from, char to)
	{
		std::replace(this->begin(), this->end(), from, to);
		return *this;
	}

	String& String::remove(char what)
	{
		this->erase(std::remove(this->begin(), this->end(), what),this->end());
		return *this;
	}

	String& String::ensureLastChar(char end)
	{
		if ( !this->hasSuffix(end) ) this->append(1, end);
		return *this;
	}

	String& String::removeWhitespaces()
	{
		Iterator end = this->end();
		end  = std::remove(this->begin(), end, ' ');
		end  = std::remove(this->begin(), end, '\t');
		end  = std::remove(this->begin(), end, '\n');
		end  = std::remove(this->begin(), end, '\r');
				
		this->erase(end,this->end());
		return *this;
	}

	String String::operator+ (int i) const
	{
		stringstream s;
		s << *this << i;
		return s.str();
	}

	String String::operator+ (unsigned int i) const
	{
		stringstream s;
		s << *this << i;
		return s.str();
	}

	String String::operator+ (short int i) const
	{
		stringstream s;
		s << *this << i;
		return s.str();
	}

	String String::operator+ (short unsigned int i) const
	{
		stringstream s;
		s << *this << i;
		return s.str();
	}

	String String::operator+ (long int i) const
	{
		stringstream s;
		s << *this << i;
		return s.str();
	}

	String String::operator+ (long unsigned int i) const
	{
		stringstream s;
		s << *this << i;
		return s.str();
	}

	String String::operator+ (long long unsigned int i) const
	{
		stringstream s;
		s << *this << i;
		return s.str();
	}

	String String::operator+ (float f) const
	{
		stringstream s;
		s.precision(7);
		s << *this << f;
		return s.str();
	}

	String String::operator+ (double d) const
	{
		stringstream s;
		s.precision(10);
		s << *this << d;
		return s.str();
	}

	String String::operator+ (long double d) const
	{
		stringstream s;
		s.precision(16);
		s << *this << d;
		return s.str();
	}

	String String::operator+ (char c) const
	{
		String tmp(*this);
		tmp.push_back(c);
		return tmp;
	}

	String String::operator+ (const char* s) const
	{
		String tmp(*this);
		tmp.append(s);
		return tmp;
	}

		String String::operator+ (const String& s) const
	{
		String tmp(*this);
		tmp.insert(tmp.end(),s.begin(),s.end());
		return tmp;
	}
		
		String String::operator+ (const std::string& s) const
	{
		String tmp(*this);
		tmp.insert(tmp.end(),s.begin(),s.end());
		return tmp;
	}

} // namespace OpenMS


