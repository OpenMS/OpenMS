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

#ifndef OPENMS_FORMAT_HANDLERS_SCHEMAHANDLER_H
#define OPENMS_FORMAT_HANDLERS_SCHEMAHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <xercesc/sax2/Attributes.hpp>

#include <map>
#include <iostream>
#include <stack>

namespace OpenMS
{
	class MetaInfoInterface;
	
	namespace Internal
	{

  /**
  	@brief abstract XML file handler supporting the use of different schemata for a format

		Works only with schemata defined in XMLSchemes.h.
		The class uses indices (or better enumeration values defined in derived classes)
		to access strings for tags or attributes of a XML-file
		(e.g. SPECTRUM instead of \<spectrum\> or \<Spectrum\>).
		This makes the implementation independent	from the underlying xml scheme.
		str2enum_() delivers the enum-value for a given string whereas enum2str_() returns
		the string for a given enum-value.
  */
  class SchemaHandler
		: public XMLHandler
  {
    public:
      /**@name Constructors and destructor */
      //@{
			
			SchemaHandler(const String& filename);
      /**
      	@brief constructor used to initialize all vectors to non-trivial sizes

				@param tag_num number of tags
				@param map_num number of maps
				@param filename the file name to handle
			*/
      SchemaHandler(Size tag_num, Size map_num,const String& filename);
      ///
      virtual ~SchemaHandler();
      //@}

			/// Finalizes members after handling a tag. Call this in your endElement() reimplementation.
			UnsignedInt leaveTag(const XMLCh* const qname);
			
			/// Sets up members for handling the current tag. Call this in your startElement() reimplementation.
			/// @returns a numerical value representing the tag.
			UnsignedInt enterTag(const XMLCh* const qname, const xercesc::Attributes& attributes);

			/// Writes the contents to a stream
			virtual void writeTo(std::ostream& os) = 0;

    protected:
    	/// Default construtctor not implemented => protected
    	SchemaHandler();
			
			std::stack<bool> skip_tag_;
			
			/// is parser currently in tag with given index?
			std::vector<bool> is_parser_in_tag_;

			/// Assoziate enumeration values with strings
			typedef std::map<std::string,int> String2EnumMap;

			/// Assoziate strings with enumeration values
			typedef std::vector<String> Enum2StringMap;

			///	vector of String2Enum-maps to map strings to enum-values
			std::vector<String2EnumMap> str2enum_array_;

			///	vector of Enum2String-maps to map an enum-value to a string
			std::vector<Enum2StringMap> enum2str_array_;

			/// index of schema from XMLSchemes.h used for Handler
			UnsignedInt schema_;

			/// pointer to attributes of current tag
			const xercesc::Attributes* atts_;
			
			/// skip parsing of the current tag
			void skipTag_();

			/// Find the enum-value that corresponds to the string @p value in map with index @p index
			UnsignedInt str2enum_(UnsignedInt index, const String& value, const char* message="");

			/** @brief Find the string that corresponds to the enum-value @p value
					in map with index @p index

					Just for convenience, in consistency with str2enum_().
			*/
			const String& enum2str_(UnsignedInt index, UnsignedInt value);

			/// Fill all str2enum-maps with strings from schema @p schema
			void fillMaps_(const String* schema);

			/// Fill particular map @p str2enum with given string array @p enum2str
			void fillMap_(String2EnumMap& str2enum, const Enum2StringMap& enum2str);

			/// Add name, value and description to a given MetaInfo object
			void setAddInfo_(	MetaInfoInterface& info, const String& name, const String& value, const String& description);
	
			/**  
				@brief write cvsParamType element containing floats to stream
				
				@p value float value
				@p acc acquisition number defined by ontology
				@p name term defined by ontology
				@p indent number of tabs used in front of tag
				
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value="@p value"/&gt;
			*/
			void writeCVS_(std::ostream& os, float value, const String& acc, const String& name, int indent=4);
	
			/**  
				@brief write cvsParamType element containing strings to stream
				
				@p value string value
				@p acc acquisition number defined by ontology
				@p name term defined by ontology
				@p indent number of tabs used in front of tag
				
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value="@p value"/&gt;
			*/
			void writeCVS_(std::ostream& os, const String& value, const String& acc, const String& name, int indent=4);
	
			/**  
				@brief write cvsParamType element containing enum-value to stream
	
				@p map index
				@p value enumeration value
				@p acc acquisition number defined by ontology
				@p name term defined by ontology
				@p indent number of tabs used in front of tag
				
				Example:
				&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value=""/&gt;
			*/
			void writeCVS_(std::ostream& os, int value, int map, const String& acc, const String& name, int indent=4);
	
			/**  
				@brief write multiple userParam elements containing MetaInfo to stream
				
				@p meta interface to access all meta info
				@p indent number of tabs used in front of tag
				
				Example:
				&lt;userParam name="??" value="??"/&gt;
			*/
			void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, int indent=4);
	
			/// check if value of attribute equals the required value, otherwise throw error
			void checkAttribute_(UnsignedInt attribute, const String& required, const String& required_alt="");
	
			/// return value of attribute as String
			String getAttributeAsString_(UnsignedInt attribute);
			
			
			void setMaps_(UnsignedInt tagmap, UnsignedInt attmap);
			
		private:
			UnsignedInt tag_map_, att_map_;
  	
  	};

	} // namespace Internal
} // namespace OpenMS

#endif
