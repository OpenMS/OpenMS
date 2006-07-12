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

#include <OpenMS/CONCEPT/Types.h>

#include <qxml.h>

namespace OpenMS
{
	namespace Internal
	{
	/**
		@brief Base class for XML handlers.
		
		This class extends the QXmlDefaultHandler by some functionality for the handling of errors.
	*/
  class XMLHandler
  	: public QXmlDefaultHandler
  {
    public:
    	/// Default constructor
      XMLHandler();
			/// Destructor
      virtual ~XMLHandler();

			/// Fatal error handler. Throws a ParseError exception
      bool fatalError(const QXmlParseException& exception);
			/// Error handler. Currently always returns false, so the parsing stops
      bool error(const QXmlParseException& exception);
			/// Warning handler. Stopts parsing depending on abort_on_warning_
      bool warning(const QXmlParseException& exception);
			
			/// Parsing method for character data
		  virtual bool characters( const QString & chars );
			/// Parsing method for opening tags
      virtual bool startElement(const QString & uri, const QString & local_name, 
																const QString & qname, const QXmlAttributes & attributes );
			/// Parsing method for closing tags
      virtual bool endElement( const QString & uri, const QString & local_name,
															 const QString & qname ); 
			
			/// Sets whether warnings are fatal
			void abortOnWarning(bool do_abort);
			/// Returns whether warnings are fatal
			bool abortOnWarning();
			
			/// Returns the last error description
  		QString errorString();
			

  	protected:
			/// Error message of the last error
			QString error_message_;
			
			/// File name
			QString file_;
			
			/// If parsing is stopped when a warning is issued
			bool abort_on_warning_;
			
			/// Conversion of a QString to an integer value
			inline SignedInt asSignedInt_(const QString& in)
			{
				bool ok = true;
				SignedInt res = in.toInt(&ok);
				if (!ok)
				{
					error(QXmlParseException(QString("SignedInt conversion error of \"%1\" parsed by %2 ").arg(in).arg(file_)));
				}
				return res;
			}
			/// Conversion of a QString to an unsigned integer value
			inline UnsignedInt asUnsignedInt_(const QString& in)
			{
				bool ok = true;
				UnsignedInt res = in.toUInt(&ok);
				if (!ok)
				{
					error(QXmlParseException(QString("UnsignedInt conversion error of \"%1\" parsed by %2 ").arg(in).arg(file_)));
				}
				return res;
			}
			/// Conversion of a QString to a double value
	 		inline double asDouble_(const QString& in)
			{
				bool ok = true;
				double res = in.toDouble(&ok);
				if (!ok)
				{
					error(QXmlParseException(QString("Double conversion error of \"%1\" parsed by %2 ").arg(in).arg(file_)));
				}
				return res;
			}
	
			/// Conversion of a QString to a float value
	 		inline float asFloat_(const QString& in)
			{
				bool ok = true;
				double res = in.toFloat(&ok);
				if (!ok)
				{
					error(QXmlParseException(QString("Float conversion error of \"%1\" parsed by %2 ").arg(in).arg(file_)));
				}
				return res;
			}
			
			/// Conversion of a QString to a bool value
	 		inline bool asBool_(const QString& in)
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
					error(QXmlParseException(QString("Boolean conversion error of \"%1\" parsed by %2 ").arg(in).arg(file_)));
				}
				return false;
			}

	};

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
