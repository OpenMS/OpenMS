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

#ifndef OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_XMLHANDLER_H

#include <iostream>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <xercesc/sax2/DefaultHandler.hpp>

namespace OpenMS
{
	namespace Internal
	{
	/**
		@brief Base class for XML handlers.
		
		@todo take advantage of XERCES string handling: XMLCh (Marc)
		@todo develop concept for optional attributes (Marc)
	*/
  class XMLHandler
  	: public xercesc::DefaultHandler
  {
    public:
    	/// Default constructor
      XMLHandler();
			/// Destructor
      virtual ~XMLHandler();

			/// Fatal error handler. Throws a ParseError exception
      void fatalError(const xercesc::SAXParseException& exception);
			/// Error handler. Currently always returns false, so the parsing stops
      void error(const xercesc::SAXParseException& exception);
			/// Warning handler.
      void warning(const xercesc::SAXParseException& exception);
			
			/// Parsing method for character data
		  virtual void characters(const XMLCh* const chars, const unsigned int length);
			/// Parsing method for opening tags
      virtual void startElement(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const xercesc::Attributes& attrs);
			/// Parsing method for closing tags
      virtual void endElement( const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname);
			
			/// Returns the last error description
  		String errorString();

  	protected:
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
					const xercesc::Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("SignedInt conversion error of \"") + in + "\" parsed by " + file_;
					error(xercesc::SAXParseException(xercesc::XMLString::transcode(message.c_str()), *loc ));
				}
				return res;
			}
			/// Conversion of a String to an unsigned integer value
			inline UnsignedInt asUnsignedInt_(const String& in)
			{
				UnsignedInt res;
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
					const xercesc::Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("UnsignedInt conversion error of \"") + in + "\" parsed by " + file_;
					error(xercesc::SAXParseException(xercesc::XMLString::transcode(message.c_str()), *loc ));
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
					const xercesc::Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("Double conversion error of \"") + in + "\" parsed by " + file_;
					error(xercesc::SAXParseException(xercesc::XMLString::transcode(message.c_str()), *loc ));
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
					const xercesc::Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("Float conversion error of \"") + in + "\" parsed by " + file_;
					error(xercesc::SAXParseException(xercesc::XMLString::transcode(message.c_str()), *loc ));
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
					const xercesc::Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("Boolean conversion error of \"") + in + "\" parsed by " + file_;
					error(xercesc::SAXParseException(xercesc::XMLString::transcode(message.c_str()), *loc ));
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
						tmp.replace('T', ' ');
						tmp = tmp.substr(0,19);
						res.set(tmp);
					}
					catch(Exception::ParseError err)
					{
						const xercesc::Locator* loc = 0;
						setDocumentLocator(loc);
						String message = String("DateTime conversion error of \"") + in + "\" parsed by " + file_;
						error(xercesc::SAXParseException(xercesc::XMLString::transcode(message.c_str()), *loc ));
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
			}

	};

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
