// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax/Locator.hpp>
#include <xercesc/sax2/Attributes.hpp>

#include <algorithm>

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
	  		: public Exception::BaseException
	  	{
	  		public:
		  		EndParsingSoftly(const char* file, int line, const char* function) 
		  			:Exception::BaseException(file,line,function)
		  		{
		  		}
	  	};

    	/// Default constructor
      XMLHandler(const String& filename, const String& version);
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
			void fatalError(const String& msg, UInt line=0, UInt column=0) const;
			/// Error handler for recoverable errors.
			void error(const String& msg, UInt line=0, UInt column=0) const;
			/// Warning handler.
			void warning(const String& msg, UInt line=0, UInt column=0) const;
			
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
			mutable String error_message_;
			
			/// File name
			String file_;
			
			/// Schema version
			String version_;
			
			/// Helper class for string conversion
			StringManager sm_;
			
			/**
				@brief Stack of open XML tags
			
				This member is used only in those XML parsers that need this information.
			*/
			std::vector<String> open_tags_;
			
			/// Returns if two xerces strings are equal
			inline bool equal_(const XMLCh* a, const XMLCh* b)
			{
				return xercesc::XMLString::compareString(a,b)==0;
			}

			///@name cvParam and userParam handling methods (for mzData and FeatureXML)
			//@{
			/**  
				@brief write cvParam containing strings to stream
				
				@p value string value
				@p acc accession number defined by ontology
				@p name term defined by ontology
				@p indent number of tabs used in front of tag
				
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value="@p value"/&gt;
			*/
			inline void writeCVS_(std::ostream& os, DoubleReal value, const String& acc, const String& name, int indent=4) const
			{
				if (value!=0.0)
				{
					os << String(indent,'\t') << "<cvParam cvLabel=\"psi\" accession=\"PSI:" << acc << "\" name=\"" << name << "\" value=\"" << value << "\"/>\n";
				}
			}
			/**  
				@brief write cvParam containing strings to stream
				
				@p value string value
				@p acc accession number defined by ontology
				@p name term defined by ontology
				@p indent number of tabs used in front of tag
				
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value="@p value"/&gt;
			*/
			inline void writeCVS_(std::ostream& os, const String& value, const String& acc, const String& name, int indent=4) const
			{
				if (value!="")
				{
					os << String(indent,'\t') << "<cvParam cvLabel=\"psi\" accession=\"PSI:" << acc << "\" name=\"" << name << "\" value=\"" << value << "\"/>\n";
				}
			}

			/**  
				@brief write cvParam element to stream

				@p os Output stream
				@p value enumeration value	
				@p map index if the terms in cv_terms_
				@p acc accession number defined by ontology
				@p name term defined by ontology
				@p indent number of tabs used in front of tag
				
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value=""/&gt;
			*/
			void writeCVS_(std::ostream& os, int value, int map, const String& acc, const String& name, int indent=4);

			///Writing the MetaInfo as UserParam to the file
			inline void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, int indent=4)
			{
				std::vector<String> keys;
				meta.getKeys(keys);
				for (std::vector<String>::const_iterator it = keys.begin(); it!=keys.end(); ++it)
				{
					if ( (*it)[0] != '#')  // internally used meta info start with '#'
					{
						os << String(indent,'\t') << "<userParam name=\"" << *it << "\" value=\"" << meta.getMetaValue(*it) << "\"/>\n";
					}
				}
			}
			//@}

			///@name General MetaInfo handling (for IdXML, FeatureXML, FeaturePairsXML)
			//@{
			///Writes the content of MetaInfoInterface to the file
			void writeUserParam_(const String& tag_name, std::ostream& os, const MetaInfoInterface& meta, UInt indent) const;
			//@}

			///@name controlled vocabulary handling methods 
			//@{
			/// Array of CV term lists (one sublist denotes one term and it's children)
			std::vector< std::vector<String> > cv_terms_;
			
			/// Converts @p term to the index of the term in the cv_terms_ entry @p section
			inline UInt cvStringToEnum_(UInt section, const String& term, const char* message)
			{
				OPENMS_PRECONDITION(section<cv_terms_.size(),"cvStringToEnum_: Index overflow (secion number too large)");
				//std::cout << "looking up key \"" << value << "\" in map nr. " << index << "..." << std::endl;
		
				std::vector<String>::const_iterator it = std::find(cv_terms_[section].begin(), cv_terms_[section].end(), term);
				if (it == cv_terms_[section].end())
				{
					std::cout << "Warning: Unexpected CV entry '" << message << "'='" << term << "' parsed in " << file_ << std::endl;
				}
				else
				{
					return  (it - cv_terms_[section].begin());
				}
				
				return 0;
			}

			//@}
			
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
			/// Conversion of a Xerces string to an integer value
			inline Int asInt_(const XMLCh* in)
			{
				return xercesc::XMLString::parseInt(in);
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
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val==0) fatalError(String("Required attribute '") + name + "' not present!");
				return sm_.convert(val);
			}
			/// Converts an attribute to a Int
			inline Int attributeAsInt_(const xercesc::Attributes& a, const char* name) const
			{
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val==0) fatalError(String("Required attribute '") + name + "' not present!");
				return xercesc::XMLString::parseInt(val);
			}
			/// Converts an attribute to a DoubleReal
			inline DoubleReal attributeAsDouble_(const xercesc::Attributes& a, const char* name) const
			{
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val==0) fatalError(String("Required attribute '") + name + "' not present!");
				return atof(sm_.convert(val));
			}
			/**
				@brief Assigns the attribute content to the String @a value if the attribute is present
				
				@return if the attribure was present
			*/
			inline bool optionalAttributeAsString_(String& value, const xercesc::Attributes& a, const char* name) const
			{
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val!=0)
				{
					char* tmp2 = sm_.convert(val);
					if (String(tmp2) != "")
					{
						value = tmp2;
						return true;
					}
				}
				return false;
			}
			/**
				@brief Assigns the attribute content to the Int @a value if the attribute is present
				
				@return if the attribure was present
			*/
			inline bool optionalAttributeAsInt_(Int& value, const xercesc::Attributes& a, const char* name) const
			{
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val!=0)
				{
					value = xercesc::XMLString::parseInt(val);
					return true;
				}
				return false;
			}
			/**
				@brief Assigns the attribute content to the UInt @a value if the attribute is present
				
				@return if the attribure was present
			*/
			inline bool optionalAttributeAsUInt_(UInt& value, const xercesc::Attributes& a, const char* name) const
			{
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val!=0)
				{
					value = xercesc::XMLString::parseInt(val);
					return true;
				}
				return false;
			}
			/**
				@brief Assigns the attribute content to the DoubleReal @a value if the attribute is present
				
				@return if the attribure was present
			*/
			inline bool optionalAttributeAsDouble_(DoubleReal& value, const xercesc::Attributes& a, const char* name) const
			{
				const XMLCh* val = a.getValue(sm_.convert(name));
				if (val!=0)
				{
					value = atof(sm_.convert(val));
					return true;
				}
				return false;
			}
			/// Converts an attribute to a String
			inline char* attributeAsString_(const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val==0) fatalError(String("Required attribute '") + sm_.convert(name) + "' not present!");
				return sm_.convert(val);
			}
			/// Converts an attribute to a Int
			inline Int attributeAsInt_(const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val==0) fatalError(String("Required attribute '") + sm_.convert(name) + "' not present!");
				return xercesc::XMLString::parseInt(val);
			}
			/// Converts an attribute to a DoubleReal
			inline DoubleReal attributeAsDouble_(const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val==0) fatalError(String("Required attribute '") + sm_.convert(name) + "' not present!");
				return atof(sm_.convert(val));
			}
			/// Assigns the attribute content to the String @a value if the attribute is present
			inline bool optionalAttributeAsString_(String& value, const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val!=0)
				{
					char* tmp2 = sm_.convert(val);
					if (String(tmp2) != "")
					{
						value = tmp2;
						return true;
					}
				}
				return false;
			}
			/// Assigns the attribute content to the Int @a value if the attribute is present
			inline bool optionalAttributeAsInt_(Int& value, const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val!=0)
				{
					value = xercesc::XMLString::parseInt(val);
					return true;
				}
				return false;
			}
			/// Assigns the attribute content to the UInt @a value if the attribute is present
			inline bool optionalAttributeAsUInt_(UInt& value, const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val!=0)
				{
					value = xercesc::XMLString::parseInt(val);
					return true;
				}
				return false;
			}
			/// Assigns the attribute content to the DoubleReal @a value if the attribute is present
			inline bool optionalAttributeAsDouble_(DoubleReal& value, const xercesc::Attributes& a, const XMLCh* name) const
			{
				const XMLCh* val = a.getValue(name);
				if (val!=0)
				{
					value = atof(sm_.convert(val));
					return true;
				}
				return false;
			}
			//@}
		
		private:
			/// Not implemented
			XMLHandler();
	};
	
	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
