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

namespace OpenMS
{
	namespace Internal
	{
	/**
		@brief Base class for XML handlers.
		
		@todo Add external CVs, integrate OLS (Marc)
		
		@todo throw out old error handlers? (Thomas K.)
	*/
  class XMLHandler
  	: public xercesc::DefaultHandler
  {
    public:
    	/// Default constructor
      XMLHandler(const String& filename);
			/// Destructor
      virtual ~XMLHandler();

			// the following three methods are reimplemented from Xerces' DefaultHandler
			// class. They forward the error message to our own error handlers below.
			/// Fatal error handler. Throws a ParseError exception
      void fatalError(const xercesc::SAXParseException& exception);
			/// Error handler. Currently always returns false, so the parsing stops
      void error(const xercesc::SAXParseException& exception);
			/// Warning handler.
      void warning(const xercesc::SAXParseException& exception);
			
			
			/// Fatal error handler. Throws a ParseError exception
			void fatalError(const String& msg);
			/// Error handler.
			void error(const String& msg);
			/// Warning handler.
			void warning(const String& msg);
			
			/// Parsing method for character data
		  virtual void characters(const XMLCh* const chars, const unsigned int length);
			/// Parsing method for opening tags
      virtual void startElement(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const xercesc::Attributes& attrs);
			/// Parsing method for closing tags
      virtual void endElement( const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname);
			
			/// Returns the last error description
  		String errorString();

  	protected:
  		XMLHandler(); /// Not implemented => protected
			
			/// Error message of the last error
			String error_message_;
			
			/// File name
			String file_;
			
			/// Conversion of a String to an integer value
			inline SignedInt asSignedInt_(const String& in)
			{
				SignedInt res = 0;
				try
				{
					res = in.toInt();
				}
				catch (Exception::ConversionError)
				{
					error(String("SignedInt conversion error of \"") + in + "\"");
				}
				return res;
			}
			/// Conversion of a String to an unsigned integer value
			inline UnsignedInt asUnsignedInt_(const String& in)
			{
				UnsignedInt res = 0;
				try
				{
					SignedInt tmp = in.toInt();
					if (tmp<0)
					{
						Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"");
					}
					res = UnsignedInt(tmp);
				}
				catch (Exception::ConversionError)
				{
					error(String("UnsignedInt conversion error of \"") + in + "\"");
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

			/// Appends the location of the @p exception to the @p message (if available) 
	 		inline void appendLocation_(const xercesc::SAXParseException& exception, String& message)
			{
				if (exception.getLineNumber()!=-1)
				{
					message = message + " at line " + String(exception.getLineNumber());
				}
				if (exception.getColumnNumber()!=-1)
				{
					message = message + " at column " + String(exception.getColumnNumber());
				}
				message = message + " in file " + file_;
			}
			
			inline void appendLocation_(const xercesc::Locator* loc, String& message)
			{
				if (loc)
				{
				  if (loc->getLineNumber() != -1) message += " at line " + String(loc->getLineNumber());
				  if (loc->getColumnNumber() != -1) message += " at column " + String(loc->getColumnNumber());
				}
				message += " in file " + file_;
			}

	};

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
