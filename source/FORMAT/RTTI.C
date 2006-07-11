// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FORMAT/RTTI.h>

#include <typeinfo>
#include <ctype.h>
#include <iostream>

// Nasty hacks to demangle the stupid name mangling schemes 
// of diverse compilers.

// GNU g++:
// Starting V3.0 we use __cxa_demangle to demangle the names,
// which is declared in <cxxabi.h>.
#ifdef OPENMS_COMPILER_GXX
#	if (OPENMS_COMPILER_VERSION_MAJOR > 2)
		#include <cxxabi.h>
#	endif
#endif

#ifdef OPENMS_COMPILER_INTEL
	// Declare the __cxa_demangle method for Intel's C++ compiler.
	// Intel does not provide the cxxabi.h header G++ provides, so 
	// this hack is somewhat rough.
	namespace abi
	{
		extern "C" char* __cxa_demangle(const char*, char*, unsigned int*, int*);
	}
#endif

namespace OpenMS 
{

	string streamClassName(const std::type_info& t)
	{
#if (defined(OPENMS_COMPILER_GXX) || defined(OPENMS_COMPILER_INTEL))
    #if (OPENMS_COMPILER_VERSION_MAJOR < 3)
			string s(t.name());
      s = GNUDemangling::demangle(s);
    #else
			char buf[OPENMS_MAX_LINE_LENGTH];
			std::size_t length = OPENMS_MAX_LINE_LENGTH - 1;
			int status = 0;
      string s("_Z");
      s += t.name();
      //std::cout <<"Demangle: "<< s << " length: "<< length<< " status: "<< status << std::endl;
      char* name = abi::__cxa_demangle(s.c_str(), buf, &length, &status);
      if (name != 0)
      {
        s = name;
			}
			//std::cout <<"End Demangle!" << std::endl;
    #endif                                                                                                                                   
#else
		string s(t.name());
		#ifdef OPENMS_COMPILER_MSVC
			// MSVC prefixes all class names with "class " -- delete it!
			while (s.find("class ") != string::npos) 
				s.erase(s.find("class "), 6);
		#endif
		
		
#endif


		for (unsigned int i = 0; i < s.size(); i++)
		{
			if (s[i] == ' ')
			{
				s[i] = '_';
			}
		}

		if (string(s, 0, 6) == "const_")
		{
			s.erase(0, 6);
		}

		return s;
	}
 
#ifdef OPENMS_COMPILER_GXX
#	if (OPENMS_COMPILER_VERSION_MAJOR < 3)

	namespace GNUDemangling 
	{

		string decode_mangling(string& s)
		{
			string tmp;
			int i,len;

			if (s.size() == 0)
				return "";

			if (!isdigit(s[0]))
			{ // decode GNU shortcuts for built-in types
				char c = s[0];
				s.erase(0, 1);
				switch (c)
				{
					case 'Q': // start of class name
						len = atoi(string(s,1,1).c_str());
						s.erase(0, 1);
						for (i = 0; i < len; i++)
						{
							tmp.append(decode_mangling(s));
							tmp.append("::");
						}
						tmp.erase(tmp.end() - 2, tmp.end());
						break;
					case 'Z': // template parameter
						return decode_mangling(s);
						break;

					case 'i':
						return "int";
						break;

					case 'l':
						return "long";
						break;

					case 's':
						return "short";
						break;

					case 'c':
						return "char";
						break;

					case 'x':
						return "long long";
						break;

					case 'f':
						return "float";
						break;

					case 'd':
						return "double";
						break;

					case 'b':
						return "bool";
						break;

					case 'w':
						return "wchar_t";
						break;

					case 'U': // unsigned variants
						tmp = "unsigned ";
						tmp.append(decode_mangling(s));
						break;

					case 'C': // const
						tmp = "const ";
						tmp.append(decode_mangling(s));
						break;

					case 'P': // pointer
						tmp = decode_mangling(s);
						tmp.append("*");
						break;

					case 'R': // reference
						tmp = decode_mangling(s);
						tmp.append("&");
						break;

					case 't':
						tmp = decode_mangling(s);
						tmp.append("<");
					
						len = atoi(string(1, s[0]).c_str());
						s.erase(0,1);
						for (i = 0; i < len; i++)
						{
							tmp.append(decode_mangling(s));
							tmp.append(",");
						}

						// remove last ','
						tmp.erase(tmp.end() - 1, tmp.end());
						tmp.append(">");
						break;

					default:
						tmp = "?";
				}
				return tmp;
			} 
			else 
			{
				i = s.find_first_not_of("0123456789");
				len = atol(string(s, 0, i).c_str());
				if (len == 0)
				{
					s.erase(0,1);
						
					if (s.size() > 0)
					{
						return decode_mangling(s);
					}
					else 
					{
						return "";
					}
				}
				else 
				{
					string h(s, i, len);
					s.erase(0, i + len);
				
					return h;
				}
			}
		}

		string demangle(string s)
		{
			string tmp = decode_mangling(s);

			while (tmp[tmp.size() - 1] == ':')
			{
				tmp.erase(tmp.end() - 1, tmp.end());
			}

			while (tmp[0] == ':')
			{
				tmp.erase(0, 1);
			}

			return tmp;
		}

	} // namespace GNUDemangling 

#	endif // (OPENMS_COMPILER_VERSION_MAJOR < 3)
#endif // OPENMS_COMPILER_GXX

} // namespace OpenMS

