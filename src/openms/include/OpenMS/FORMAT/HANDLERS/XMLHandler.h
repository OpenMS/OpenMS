// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Macros.h>

 
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/FORMAT/HANDLERS/StringManager.h>
#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>

#include <iosfwd>
#include <string>
#include <string_view>


namespace OpenMS
{
  class ControlledVocabulary;
  class CVTerm;
  class MetaInfoInterface;
  class ProteinIdentification;

  namespace Internal
  {

    // XMLCh string handling lives in StringManager (Xerces-free public API);
    // shared_xerces_ptr / unique_xerces_ptr / CONST_XMLCH are declared there.

    /**
        @brief Base class for XML handlers.
    */
    class OPENMS_DLLAPI XMLHandler
    {
public:

      /// Exception that is thrown if the parsing is ended by some event (e.g. if only a prefix of the XML file is needed).
      class OPENMS_DLLAPI EndParsingSoftly :
        public Exception::BaseException
      {
      public:
        EndParsingSoftly(const char * file, int line, const char * function) :
          Exception::BaseException(file, line, function)
        {
        }

      };

      ///Action to set the current mode (for error messages)
      enum ActionMode
      {
        LOAD,               ///< Loading a file
        STORE               ///< Storing a file
      };

      enum LOADDETAIL 
      {  
        LD_ALLDATA,       // default; load all data
        LD_RAWCOUNTS,     // only count the total number of spectra and chromatograms (usually very fast)
        LD_COUNTS_WITHOPTIONS // count the number of spectra, while respecting PeakFileOptions (msLevel and RTRange) and chromatograms (fast)
      };


      /// Default constructor
      XMLHandler(const std::string & filename, const std::string & version);
      /// Destructor
      virtual ~XMLHandler();

      /// Release internal memory used for parsing (call
      void reset();

      /// Fatal error handler. Throws a ParseError exception
      void fatalError(ActionMode mode, const std::string & msg, UInt line = 0, UInt column = 0) const;
      /// Error handler for recoverable errors.
      void error(ActionMode mode, const std::string & msg, UInt line = 0, UInt column = 0) const;
      /// Warning handler.
      void warning(ActionMode mode, const std::string & msg, UInt line = 0, UInt column = 0) const;

      /**
          @name Native (Xerces-free) SAX callbacks

          Subclasses override these. @c qname / @c chars are UTF-16 code units (as
          delivered by the parser); attributes are exposed through the Xerces-free
          @ref XMLAttributes view. The Xerces SAX parser is bridged to these
          callbacks by the internal SAX2HandlerAdapter.
      */
      //@{
      /// Parsing method for opening tags
      virtual void onStartElement(const char16_t * qname, const XMLAttributes & attributes);
      /// Parsing method for closing tags
      virtual void onEndElement(const char16_t * qname);
      /// Parsing method for character data
      virtual void onCharacters(const char16_t * chars, Size length);
      //@}

      /// Writes the contents to a stream.
      virtual void writeTo(std::ostream & /*os*/);

      /// handler which support partial loading, implement this method
      virtual LOADDETAIL getLoadDetail() const;

      /// handler which support partial loading, implement this method
      virtual void setLoadDetail(const LOADDETAIL d);

      /**
        @brief Escapes a string and returns the escaped string

        Some characters must be escaped which are allowed in user params. E.g. > and & are not in XML and
        need to be escaped. Parsing those escaped strings from file again is automatically done by Xerces.
        Escaped characters are: & < > " ' 
      */
      static std::string writeXMLEscape(const std::string& to_escape)
      {
        std::string _copy = to_escape;
        // has() is cheap, so check before calling substitute(), since substitute() will usually happen rarely
        if (StringUtils::has(_copy, '&')) StringUtils::substitute(_copy, "&","&amp;");
        if (StringUtils::has(_copy, '>')) StringUtils::substitute(_copy, ">","&gt;");
        if (StringUtils::has(_copy, '"')) StringUtils::substitute(_copy, "\"","&quot;");
        if (StringUtils::has(_copy, '<')) StringUtils::substitute(_copy, "<","&lt;");
        if (StringUtils::has(_copy, '\'')) StringUtils::substitute(_copy, "'","&apos;");

        return _copy;
      }

      /**
      *  @brief Convert an XSD type (e.g. 'xsd:double') to a DataValue.
      * 
      *  Not all conversions are supported yet, due to DataValue using an Int64 as the largest possible integer.
      *  Thus, if the @p value contains a large UInt64, conversion will fail.
      *  Value ranges are currently also not checked, only for XSD types which happen to match the internal representation.
      * 
      *  @param[out] type An XSD type. If the type is not supported, the returned type will be a string
      *  @param[in] value The value in sting format, e.g. "123.34"
      *  @return The Datavalue with the respective type (double, int or string)
      *  @throws Exception::ConversionError if the value does not fit into the internal representation or (for few types) exceeds the XSD specs.
      * 
      */
      static DataValue fromXSDString(const std::string& type, const std::string& value)
      {
        DataValue data_value;
        // float type
        if (type == "xsd:double" || type == "xsd:float" || type == "xsd:decimal")
        {
          data_value = DataValue(StringUtils::toDouble(value));
        }
        // <=32 bit integer types
        else if (type == "xsd:byte" ||          // 8bit signed
                 type == "xsd:int" ||           // 32bit signed
                 type == "xsd:unsignedShort" || // 16bit unsigned
                 type == "xsd:short" ||         // 16bit signed
                 type == "xsd:unsignedByte" || type == "xsd:unsignedInt")
        {
          data_value = DataValue(StringUtils::toInt32(value));
        }
        // 64 bit integer types
        else if (type == "xsd:long" || type == "xsd:unsignedLong" ||       // 64bit signed or unsigned respectively
                 type == "xsd:integer" || type == "xsd:negativeInteger" || // any 'integer' has arbitrary size... but we have to cope with 64bit for now.
                 type == "xsd:nonNegativeInteger" || type == "xsd:nonPositiveInteger" || type == "xsd:positiveInteger")
        {
          data_value = DataValue(StringUtils::toInt64(value)); // internally a signed 64-bit integer. So if someone uses 2^64-1 as value, toInt64() will raise an exception...
        }
        // everything else is treated as a string
        else
        {
          data_value = DataValue(value);
        }
        return data_value;
      }


      /**
         @brief Convert the value of a <em>\<cvParam value=.\></em> (as commonly found in PSI schemata) to the DataValue with the correct type (e.g. int) according to
                the type stored in the CV (usually PSI-MS CV), as well as set its unit.

         @param[in] cv A CV, usually the PSI-MS CV, see ControlledVocabulary::getPSIMSCV()
         @param[in] parent_tag The tag which encloses the \<cvParam\>
         @param[in] accession The accession from the 'accession' attribute of the \<cvParam\>
         @param[in] name The name from the 'name' attribute of the \<cvParam\>
         @param[in] value The value from the 'value' attribute of the \<cvParam\>
         @param[in] unit_accession The unit_accession from the 'unitAccession' attribute of the \<cvParam\>
         @return DataValue::EMPTY if a conversion error occured (e.g. if @p value could not be converted to an integer for an @p accession which requires an integer) or the DataValue upon success
      */
      DataValue cvParamToValue(const ControlledVocabulary& cv, const std::string& parent_tag, 
                               const std::string& accession, const std::string& name, const std::string& value,
                               const std::string& unit_accession) const;

      /**
         @brief Convert the value of a <em>\<cvParam value=.\></em> (as commonly found in PSI schemata) to the DataValue with the correct type (e.g. int) according to
                the type stored in the CV (usually PSI-MS CV), as well as set its unit.

         @param[in] cv A CV, usually the PSI-MS CV, see ControlledVocabulary::getPSIMSCV()
         @param[in] raw_term Represenation of the raw data (i.e. all strings) from a \<cvParam ...\> without the conversion to a specific value type
         @return DataValue::EMPTY if a conversion error occured (e.g. if @p value could not be converted to an integer for an @p accession which requires an integer) or the DataValue upon success
      */
      DataValue cvParamToValue(const ControlledVocabulary& cv, const CVTerm& raw_term) const;

      /// throws a ParseError if protIDs are not unique, i.e. PeptideIDs will be randomly assigned (bad!)
      /// Should be called before writing any ProtIDs to file
      void checkUniqueIdentifiers_(const std::vector<ProteinIdentification>& prot_ids) const;

protected:
      /// File name
      std::string file_;

      /// Schema version
      std::string version_;

      /// Helper class for string conversion
      StringManager sm_;

      /**
          @brief Stack of open XML tags

          This member is used only in those XML parsers that need this information.
      */
      std::vector<std::string> open_tags_;

      /// parse only until total number of scans and chroms have been determined from attributes
      LOADDETAIL load_detail_; 


      /// Returns if two Xerces strings are equal
      inline bool equal_(const char16_t * a, const char16_t * b) const
      {
        // Guard null pointers before constructing views: u16string_view(nullptr)
        // is UB, whereas the previous XMLString::compareString tolerated nulls
        // (both null => equal, one null => unequal), which this reproduces.
        if (a == nullptr || b == nullptr) return a == b;
        return std::u16string_view(a) == std::u16string_view(b);
      }

      ///@name General MetaInfo handling (for idXML, featureXML, consensusXML)
      //@{

      /// Writes the content of MetaInfoInterface to the file
      void writeUserParam_(const std::string & tag_name, std::ostream & os, const MetaInfoInterface & meta, UInt indent) const;

      //@}

      ///@name controlled vocabulary handling methods
      //@{

      /// Array of CV term lists (one sublist denotes one term and it's children)
      std::vector<std::vector<std::string> > cv_terms_;

      /// Converts @p term to the index of the term in the cv_terms_ entry @p section
      /// If the term is not found, @p result_on_error is returned (0 by default)
      SignedSize cvStringToEnum_(const Size section, const std::string & term, const char * message, const SignedSize result_on_error = 0);

      //@}

      ///@name std::string conversion
      //@{

      /// Conversion of a std::string to an integer value
      inline Int asInt_(const std::string & in) const
      {
        Int res = 0;
        try
        {
          res = StringUtils::toInt32(in);
        }
        catch (Exception::ConversionError&)
        {
          error(LOAD,std::string("Int conversion error of \"") + in + "\"");
        }
        return res;
      }

      /// Conversion of a Xerces string to an integer value
      inline Int asInt_(const char16_t * in) const
      {
        return sm_.parseInt(in);
      }

      /// Conversion of a std::string to an unsigned integer value
      inline UInt asUInt_(const std::string & in) const
      {
        UInt res = 0;
        try
        {
          Int tmp = StringUtils::toInt32(in);
          if (tmp < 0)
          {
            throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "");
          }
          res = UInt(tmp);
        }
        catch (Exception::ConversionError& )
        {
          error(LOAD,std::string("UInt conversion error of \"") + in + "\"");
        }
        return res;
      }

      /// Conversion of a std::string to a double value
      inline double asDouble_(const std::string & in) const
      {
        double res = 0.0;
        try
        {
          res = StringUtils::toDouble(in);
        }
        catch (Exception::ConversionError& )
        {
          error(LOAD,std::string("Double conversion error of \"") + in + "\"");
        }
        return res;
      }

      /// Conversion of a std::string to a float value
      inline float asFloat_(const std::string & in) const
      {
        float res = 0.0;
        try
        {
          res = StringUtils::toFloat(in);
        }
        catch (Exception::ConversionError& )
        {
          error(LOAD,std::string("Float conversion error of \"") + in + "\"");
        }
        return res;
      }

      /**
          @brief Conversion of a string to a boolean value

          'true', 'false', '1' and '0' are accepted.

          @n For all other values a parse error is produced.
      */
      inline bool asBool_(const std::string & in) const
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
          error(LOAD,std::string("Boolean conversion error of \"") + in + "\"");
        }
        return false;
      }

      /// Conversion of a xs:datetime string to a DateTime value
      inline DateTime asDateTime_(std::string date_string) const
      {
        DateTime date_time;
        if (!date_string.empty())
        {
          try
          {
            //strip away milliseconds
            StringUtils::trim(date_string);
            date_string = StringUtils::substr(date_string, 0, 19);
            date_time.set(date_string);
          }
          catch (Exception::ParseError& /*err*/ )
          {
            error(LOAD,std::string("DateTime conversion error of \"") + date_string + "\"");
          }
        }
        return date_time;
      }

      //@}

      ///@name Accessing attributes
      //@{

      /// Converts an attribute to a String
      inline std::string attributeAsString_(const XMLAttributes & a, const char * name) const
      {
        const char16_t * val = a.value(name);
        if (val == nullptr) fatalError(LOAD,std::string("Required attribute '") + name + "' not present!");
        return sm_.convert(val);
      }

      /// Converts an attribute to a Int
      inline Int attributeAsInt_(const XMLAttributes & a, const char * name) const
      {
        const char16_t * val = a.value(name);
        if (val == nullptr) fatalError(LOAD,std::string("Required attribute '") + name + "' not present!");
        return sm_.parseInt(val);
      }

      /// Converts an attribute to a double
      inline double attributeAsDouble_(const XMLAttributes & a, const char * name) const
      {
        const char16_t * val = a.value(name);
        if (val == nullptr) fatalError(LOAD,std::string("Required attribute '") + name + "' not present!");
        return StringUtils::toDouble(sm_.convert(val));
      }

      /// Converts an attribute to a DoubleList
      DoubleList attributeAsDoubleList_(const XMLAttributes & a, const char * name) const;

      /// Converts an attribute to an IntList
      IntList attributeAsIntList_(const XMLAttributes & a, const char * name) const;

      /// Converts an attribute to an StringList
      StringList attributeAsStringList_(const XMLAttributes & a, const char * name) const;

      /**
          @brief Assigns the attribute content to the String @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsString_(std::string & value, const XMLAttributes & a, const char * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = sm_.convert(val);
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the Int @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsInt_(Int & value, const XMLAttributes & a, const char * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = sm_.parseInt(val);
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the UInt @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsUInt_(UInt & value, const XMLAttributes & a, const char * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = sm_.parseInt(val);
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the double @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsDouble_(double & value, const XMLAttributes & a, const char * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value =StringUtils::toDouble(sm_.convert(val));
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the DoubleList @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsDoubleList_(DoubleList & value, const XMLAttributes & a, const char * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = attributeAsDoubleList_(a, name);
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the StringList @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsStringList_(StringList & value, const XMLAttributes & a, const char * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = attributeAsStringList_(a, name);
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the IntList @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsIntList_(IntList & value, const XMLAttributes & a, const char * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = attributeAsIntList_(a, name);
          return true;
        }
        return false;
      }

      /// Converts an attribute to a String
      inline std::string attributeAsString_(const XMLAttributes & a, const char16_t * name) const
      {
        const char16_t * val = a.value(name);
        if (val == nullptr) fatalError(LOAD,std::string("Required attribute '") + sm_.convert(name) + "' not present!");
        return sm_.convert(val);
      }

      /// Converts an attribute to a Int
      inline Int attributeAsInt_(const XMLAttributes & a, const char16_t * name) const
      {
        const char16_t * val = a.value(name);
        if (val == nullptr) fatalError(LOAD,std::string("Required attribute '") + sm_.convert(name) + "' not present!");
        return sm_.parseInt(val);
      }

      /// Converts an attribute to a double
      inline double attributeAsDouble_(const XMLAttributes & a, const char16_t * name) const
      {
        const char16_t * val = a.value(name);
        if (val == nullptr) fatalError(LOAD,std::string("Required attribute '") + sm_.convert(name) + "' not present!");
        return StringUtils::toDouble(sm_.convert(val));
      }

      /// Converts an attribute to a DoubleList
      DoubleList attributeAsDoubleList_(const XMLAttributes & a, const char16_t * name) const;

      /// Converts an attribute to a IntList
      IntList attributeAsIntList_(const XMLAttributes & a, const char16_t * name) const;

      /// Converts an attribute to a StringList
      StringList attributeAsStringList_(const XMLAttributes & a, const char16_t * name) const;

      /// Assigns the attribute content to the String @a value if the attribute is present
      inline bool optionalAttributeAsString_(std::string& value, const XMLAttributes & a, const char16_t * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = sm_.convert(val);
          return !value.empty();
        }
        return false;
      }

      /// Assigns the attribute content to the Int @a value if the attribute is present
      inline bool optionalAttributeAsInt_(Int & value, const XMLAttributes & a, const char16_t * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = sm_.parseInt(val);
          return true;
        }
        return false;
      }

      /// Assigns the attribute content to the UInt @a value if the attribute is present
      inline bool optionalAttributeAsUInt_(UInt & value, const XMLAttributes & a, const char16_t * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = sm_.parseInt(val);
          return true;
        }
        return false;
      }

      /// Assigns the attribute content to the double @a value if the attribute is present
      inline bool optionalAttributeAsDouble_(double & value, const XMLAttributes & a, const char16_t * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = StringUtils::toDouble(sm_.convert(val));
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the DoubleList @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsDoubleList_(DoubleList & value, const XMLAttributes & a, const char16_t * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = attributeAsDoubleList_(a, name);
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the IntList @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsIntList_(IntList & value, const XMLAttributes & a, const char16_t * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = attributeAsIntList_(a, name);
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the StringList @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsStringList_(StringList & value, const XMLAttributes & a, const char16_t * name) const
      {
        const char16_t * val = a.value(name);
        if (val != nullptr)
        {
          value = attributeAsStringList_(a, name);
          return true;
        }
        return false;
      }

      //@}

private:
      /// Not implemented
      XMLHandler();

      inline const std::string& expectList_(const std::string& str) const
      {
        if (!(StringUtils::hasPrefix(str, '[') && StringUtils::hasSuffix(str, ']')))
        {
          fatalError(LOAD,std::string("List argument is not a string representation of a list!"));
        }
        return str;
      }

    };

  }   // namespace Internal
} // namespace OpenMS


