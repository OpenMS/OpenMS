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
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_SCHEMAHANDLER_H
#define OPENMS_FORMAT_HANDLERS_SCHEMAHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <map>
#include <qstring.h>

namespace OpenMS
{
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
			///
			SchemaHandler();
      /** @brief constructor used to initialize all vectors to non-trivial sizes

				@p tag_num number of tags
				@p map_num number of maps
			*/
      SchemaHandler(Size tag_num, Size map_num);
      ///
      virtual ~SchemaHandler();
      //@}

    protected:
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

			/// Find the enum-value that corresponds to the string @p value in map
			/// with index @p index
			inline UnsignedInt str2enum_(UnsignedInt index, QString value, const char* message="")
			{
				String2EnumMap::const_iterator it =  str2enum_array_[index].find(value.ascii());
				if (it == str2enum_array_[index].end()){  // no enum-value for string defined
					if (useWarnings())
						warning(QXmlParseException(QString("Unhandled %3 \"%1\" parsed by %2")
																				.arg(value).arg(file_).arg(message)));
				}	else
					return it->second;
				return 0;
			}

			/** @brief Find the string that corresponds to the enum-value @p value
					in map with index @p index

					Just for convenience, in consistency with str2enum_().
			*/
			inline String enum2str_(UnsignedInt index, UnsignedInt value)
			{
				return enum2str_array_[index][value];
			}

			/// Fill all str2enum-maps with strings from schema @p schema
			void fillMaps_(const String* schema)
			{
				for (Size i= 0; i<str2enum_array_.size(); i++)
				{
					//i=0 contains scheme name -> i+1
					schema[i+1].split(';',enum2str_array_[i]);
					fillMap_(str2enum_array_[i], enum2str_array_[i]);
				}
			}

			/// Fill particular map @p str2enum with given string array @p enum2str
			inline void fillMap_(String2EnumMap& str2enum, const Enum2StringMap& enum2str)
			{
				for (Size i=0; i<enum2str.size(); i++)
				{
					str2enum[ enum2str[i] ] = i;
				}
			}

		/// Add name, value and description to a given MetaInfo object
		inline void setAddInfo(	MetaInfoInterface& info, QString name,
														QString value, String description)
		{
			info.metaRegistry().registerName(name.ascii(), description);
			info.setMetaValue(name.ascii(),value.ascii());
		}

		/**  @brief write cvsParamType element containing floats to stream */
		/**
				 @p value float value
				 @p acc acquisition number defined by ontology
				 @p name term defined by ontology
				 @p indent number of tabs used in front of tag

					Example:
					&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value="@p value"/&gt;
		*/
		inline void writeCVS_(std::ostream& os, float value, String acc, String name, int indent=4)
		{
			if (value)
				os << String(indent,'\t') << "<cvParam cvLabel=\"psi\" accession=\"PSI:"
					 << acc << "\" name=\""
					 << name << "\" value=\""
					 << value << "\"/>\n";
		}

		/**  @brief write cvsParamType element containing strings to stream */
		/**
				 @p value string value
				 @p acc acquisition number defined by ontology
				 @p name term defined by ontology
				 @p indent number of tabs used in front of tag

					Example:
					&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value="@p value"/&gt;
		*/
		inline void writeCVS_(std::ostream& os, String value, String acc, String name, int indent=4)
		{
			if (value!="")
				os << String(indent,'\t') << "<cvParam cvLabel=\"psi\" accession=\"PSI:"
					 << acc << "\" name=\""
					 << name << "\" value=\""
					 << value << "\"/>\n";
		}

		/**  @brief write cvsParamType element containing enum-value to stream */
		/**
				 @p map index
				 @p value enumeration value
				 @p acc acquisition number defined by ontology
				 @p name term defined by ontology
				 @p indent number of tabs used in front of tag

					Example:
					&lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value=""/&gt;
		*/
		inline void writeCVS_(std::ostream& os, int value, int map,
													String acc, String name, int indent=4)
		{
			writeCVS_(os, enum2str_(map,value), acc, name, indent);
		}

		/**  @brief write multiple userParam elements containing MetaInfo to stream */
		/**
				 @p meta interface to access all meta info
				 @p indent number of tabs used in front of tag

					Example:
					&lt;userParam name="??" value="??"/&gt;
		*/
		inline void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, int indent=4)
		{
			std::vector<std::string> keys;  // Vector to hold keys to meta info
			meta.getKeys(keys);

			for (std::vector<std::string>::const_iterator it = keys.begin(); it!=keys.end(); ++it)
				if ( (*it)[0] != '#')  // internally used meta info start with '#'
					os << String(indent,'\t') << "<userParam name=\""
						 << *it << "\" value=\""
						 << meta.getMetaValue(*it) << "\"/>\n";
		}

  };

	} // namespace Internal
} // namespace OpenMS

#endif
