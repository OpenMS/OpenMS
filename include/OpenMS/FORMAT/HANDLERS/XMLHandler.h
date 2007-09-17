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

#ifndef OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_XMLHANDLER_H

#include <iostream>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax/Locator.hpp>
#include <xercesc/sax2/Attributes.hpp>

namespace OpenMS
{
	namespace Internal
	{
	
		/// Helper class for XML parsing that handles the memory management for conversions of Xerces strings
		class StringManager 
		{
			public:
				/// Constructor
				StringManager();
			
				/// Destructor. Destroys the strings in the various lists
				~StringManager();
						
				/// Frees memory of all owned strings
				void clear();
			
				/// Transcode the supplied C string to XMLCh* and take ownership of the XMLCh*
				XMLCh* convert(const char* str) const;
			
				/// Transcode the supplied C++ string to XMLCh* and take ownership of the XMLCh*			
				XMLCh* convert(const std::string& str) const;

				/// Transcode the supplied OpenMS string to XMLCh* and take ownership of the XMLCh*			
				XMLCh* convert(const String& str) const;
			
				/// Transcode the supplied XMLCh* to a C string and take ownership of the C string
				char* convert(const XMLCh* str) const;
			private:
				mutable std::vector<XMLCh*> xml_strings_ ;
				mutable std::vector<char*> c_strings_ ;			
		};

	/**
		@brief Base class for XML handlers.
	*/
  class XMLHandler
  	: public xercesc::DefaultHandler
  {
    public:
	  	/// Exception that is thrown if the parsing is ended by some event (e.g. if only a prefix of the XML file is needed).
	  	class EndParsingSoftly
	  		: public Exception::Base
	  	{
	  		public:
		  		EndParsingSoftly(const char* file, int line, const char* function) throw()
		  			:Exception::Base(file,line,function)
		  		{
		  		}
	  	};

    	/// Default constructor
      XMLHandler(const String& filename);
			/// Destructor
      virtual ~XMLHandler();

			/**
				@name Reimplemented XERCES-C error handlers
				
				These methods forward the error message to our own error handlers below.
			*/
			//@{
			void fatalError(const xercesc::SAXParseException& exception);
			void error(const xercesc::SAXParseException& exception);
			void warning(const xercesc::SAXParseException& exception);
			//@}
			
			/// Fatal error handler. Throws a ParseError exception
			void fatalError(const String& msg, UInt line=0, UInt column=0);
			/// Error handler for recoverable errors.
			void error(const String& msg, UInt line=0, UInt column=0);
			/// Warning handler.
			void warning(const String& msg, UInt line=0, UInt column=0);
			
			/// Parsing method for character data
		  virtual void characters(const XMLCh* const chars, unsigned int length);
			/// Parsing method for opening tags
      virtual void startElement(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const xercesc::Attributes& attrs);
			/// Parsing method for closing tags
      virtual void endElement( const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname);

			/// Writes the contents to a stream.
			virtual void writeTo(std::ostream& /*os*/) {};
			
			/// Returns the last error description
  		String errorString();

  	protected:
			/// Error message of the last error
			String error_message_;
			
			/// File name
			String file_;
			
			/// Helper class for string conversion
			StringManager sm_;
			
			/// Returns if two xerces strings are equal
			inline bool equal(const XMLCh* a, const XMLCh* b)
			{
				return xercesc::XMLString::compareString(a,b)==0;
			}

			
			///@name String conversion
			//@{ 
			/// Conversion of a String to an integer value
			inline Int asInt_(const String& in)
			{
				Int res = 0;
				try
				{
					res = in.toInt();
				}
				catch (Exception::ConversionError)
				{
					error(String("Int conversion error of \"") + in + "\"");
				}
				return res;
			}
			/// Conversion of a String to an unsigned integer value
			inline UInt asUInt_(const String& in)
			{
				UInt res = 0;
				try
				{
					Int tmp = in.toInt();
					if (tmp<0)
					{
						Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"");
					}
					res = UInt(tmp);
				}
				catch (Exception::ConversionError)
				{
					error(String("UInt conversion error of \"") + in + "\"");
				}
				return res;
			}
			/// Conversion of a String to a double value
	 		inline double asDouble_(const String& in)
			{
				double res = 0.0;
				try
				{
					res = in.toDouble();
				}
				catch (Exception::ConversionError)
				{
					error(String("Double conversion error of \"") + in + "\"");
				}
				return res;
			}
			/// Conversion of a String to a float value
	 		inline float asFloat_(const String& in)
			{
				float res = 0.0;
				try
				{
					res = in.toFloat();
				}
				catch (Exception::ConversionError)
				{
					error(String("Float conversion error of \"") + in + "\"");
				}
				return res;
			}
			/// Conversion of a String to a bool value
	 		inline bool asBool_(const String& in)
			{
				if (in == "true" || in == "TRUE" || in == "True" || in == "1") 
				{
					return true;
				}
				else if (in == "false" || in == "FALSE" || in == "False" || in == "0")
				{
					 return false;
				}
				else 
				{
					error(String("Boolean conversion error of \"") + in + "\"");
				}
				return false;
			}
			/// Conversion of a String to a DataTime value
	 		inline DateTime asDateTime_(const String& in)
			{
				//std::cout << "IN: " << in << std::endl;
				DateTime res;
				if (in!="")
				{ 
					try
					{
						// xs::DateTime to OpenMS::DateTime
						String tmp(in);
						tmp.substitute('T', ' ');
						tmp = tmp.substr(0,19);
						res.set(tmp);
					}
					catch(Exception::ParseError err)
					{
						error(String("DateTime conversion error of \"") + in + "\"");
					}
				}
				return res;
			}
			//@}
		
			///@name Accessing attributes
			//@{		
			/// Converts an attribute to a String
			inline char* attributeAsString_(const xercesc::Attributes& a, const char* name) const
			{
				return sm_.convert(a.getValue(sm_.convert(name)));
			}
			/// Converts an attribute to a Int
			inline Int attributeAsInt_(const xercesc::Attributes& a, const char* name) const
			{
				return xercesc::XMLString::parseInt(a.getValue(sm_.convert(name)));
			}
			/// Converts an attribute to a DoubleReal
			inline DoubleReal attributeAsDouble_(const xercesc::Attributes& a, const char* name) const
			{
				return atof(sm_.convert(a.getValue(sm_.convert(name))));
			}
			/// Assigns the attribute content to the String @a value if the attribute is present
			inline void optionalAttributeAsString_(String& value, const xercesc::Attributes& a, const char* name) const
			{
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val!=0)
				{
					char* tmp2 = sm_.convert(val);
					if (String(tmp2) != "")
					{
						value = tmp2;
					}
				}
			}
			/// Assigns the attribute content to the Int @a value if the attribute is present
			inline void optionalAttributeAsInt_(Int& value, const xercesc::Attributes& a, const char* name) const
			{
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val!=0)
				{
					value = xercesc::XMLString::parseInt(val);
				}
			}
			/// Assigns the attribute content to the UInt @a value if the attribute is present
			inline void optionalAttributeAsUInt_(UInt& value, const xercesc::Attributes& a, const char* name) const
			{
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val!=0)
				{
					value = xercesc::XMLString::parseInt(val);
				}
			}
			/// Assigns the attribute content to the DoubleReal @a value if the attribute is present
			inline void optionalAttributeAsDouble_(DoubleReal& value, const xercesc::Attributes& a, const char* name) const
			{
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val!=0)
				{
					value = atof(sm_.convert(val));
				}
			}
			/// Converts an attribute to a String
			inline char* attributeAsString_(const xercesc::Attributes& a, const XMLCh* name) const
			{
				return sm_.convert(a.getValue(name));
			}
			/// Converts an attribute to a Int
			inline Int attributeAsInt_(const xercesc::Attributes& a, const XMLCh* name) const
			{
				return xercesc::XMLString::parseInt(a.getValue(name));
			}
			/// Converts an attribute to a DoubleReal
			inline DoubleReal attributeAsDouble_(const xercesc::Attributes& a, const XMLCh* name) const
			{
				return atof(sm_.convert(a.getValue(name)));
			}
			/// Assigns the attribute content to the String @a value if the attribute is present
			inline void optionalAttributeAsString_(String& value, const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val!=0)
				{
					char* tmp2 = sm_.convert(val);
					if (String(tmp2) != "")
					{
						value = tmp2;
					}
				}
			}
			/// Assigns the attribute content to the Int @a value if the attribute is present
			inline void optionalAttributeAsInt_(Int& value, const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val!=0)
				{
					value = xercesc::XMLString::parseInt(val);
				}
			}
			/// Assigns the attribute content to the UInt @a value if the attribute is present
			inline void optionalAttributeAsUInt_(UInt& value, const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val!=0)
				{
					value = xercesc::XMLString::parseInt(val);
				}
			}
			/// Assigns the attribute content to the DoubleReal @a value if the attribute is present
			inline void optionalAttributeAsDouble_(DoubleReal& value, const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val!=0)
				{
					value = atof(sm_.convert(val));
				}
			}
			//@}
		
		private:
			/// Not implemented
			XMLHandler();
	};
	
	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
