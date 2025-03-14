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

#include <OpenMS/DATASTRUCTURES/ListUtils.h> // StringList
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/util/XMLString.hpp>

#include <iosfwd>
#include <string>
#include <memory>


namespace OpenMS
{
  class ControlledVocabulary;
  class CVTerm;
  class MetaInfoInterface;
  class ProteinIdentification;

  namespace Internal
  {

    #define CONST_XMLCH(s) reinterpret_cast<const ::XMLCh*>(u ## s)

    static_assert(sizeof(::XMLCh) == sizeof(char16_t),
                  "XMLCh is not sized correctly for UTF-16.");

    //Adapted from https://www.codeproject.com/articles/99551/redux-raii-adapter-for-xerces
    //Copyright 2010 Orjan Westin
    //Under BSD license
    //========================================================================================================
    template<typename T>
    class OPENMS_DLLAPI shared_xerces_ptr
    {
      // Function to release Xerces data type with a release member function
      template<typename U>
      static void doRelease_(U* item)
      {
        // Only release this if it has no owner
        if (nullptr == item->getOwnerDocument())
          item->release();
      }

      static void doRelease_(char* item);
      static void doRelease_(XMLCh* item);

      // The actual data we're holding
      std::shared_ptr<T> item_;
    public:
      // Default constructor
      shared_xerces_ptr() = default;
      // Assignment constructor
      shared_xerces_ptr(T* item)
          : item_(item, doRelease_ )
      {}
      // Assignment of data to guard
      shared_xerces_ptr& operator=(T* item)
      {
        assign(item);
        return *this;
      }
      // Give up hold on data
      void reset()
      {
        item_.reset();
      }
      // Release currently held data, if any, to hold another
      void assign(T* item)
      {
        item_.reset(item, doRelease_ );
      }
      // Get pointer to the currently held data, if any
      T* get()
      {
        return item_.get();
      }
      const T* get() const
      {
        return item_.get();
      }
      // Return true if no data is held
      bool is_released() const
      {
        return (nullptr == item_.get());
      }
    };

    template <typename T>
    class OPENMS_DLLAPI unique_xerces_ptr
    {
    private:

      template<typename U>
      static void doRelease_(U*& item)
      {
        // Only release this if it has no parent (otherwise
        // parent will release it)
        if (nullptr == item->getOwnerDocument())
          item->release();
      }

      static void doRelease_(char*& item);
      static void doRelease_(XMLCh*& item);

      T* item_;

    public:

      // Hide copy constructor and assignment operator
      unique_xerces_ptr(const unique_xerces_ptr<T>&) = delete;
      unique_xerces_ptr& operator=(const unique_xerces_ptr<T>&) = delete;

      unique_xerces_ptr()
          : item_(nullptr)
      {}

      explicit unique_xerces_ptr(T* i)
          : item_(i)
      {}

      ~unique_xerces_ptr()
      {
        xerces_release();
      }

      unique_xerces_ptr(unique_xerces_ptr<T>&& other) noexcept
          : item_(nullptr)
      {
        this->swap(other);
      }

      void swap(unique_xerces_ptr<T>& other) noexcept
      {
        std::swap(item_, other.item_);
      }

      // Assignment of data to guard (not chainable)
      void operator=(T* i)
      {
        reassign(i);
      }

      // Release held data (i.e. delete/free it)
      void xerces_release()
      {
        if (!is_released())
        {
          // Use type-specific release mechanism
          doRelease_(item_);
          item_ = nullptr;
        }
      }

      // Give up held data (i.e. return data without releasing)
      T* yield()
      {
        T* tempItem = item_;
        item_ = nullptr;
        return tempItem;
      }

      // Release currently held data, if any, to hold another
      void assign(T* i)
      {
        xerces_release();
        item_ = i;
      }

      // Get pointer to the currently held data, if any
      T* get() const
      {
        return item_;
      }

      // Return true if no data is held
      bool is_released() const
      {
        return (nullptr == item_);
      }
    };

    //========================================================================================================

    /*
     * @brief Helper class for XML parsing that handles the conversions of Xerces strings
     *
     * It provides the convert() function which internally calls
     * XMLString::transcode and ensures that the memory is released properly
     * through XMLString::release internally. It returns a std::string or
     * std::basic_string<XMLCh> to the caller who takes ownership of the data.
     *
    */
    class OPENMS_DLLAPI StringManager
    {

      typedef std::basic_string<XMLCh> XercesString;

      // Converts from a narrow-character string to a wide-character string.
      inline static unique_xerces_ptr<XMLCh> fromNative_(const char* str)
      {
        return unique_xerces_ptr<XMLCh>(xercesc::XMLString::transcode(str));
      }

      // Converts from a narrow-character string to a wide-character string.
      inline static unique_xerces_ptr<XMLCh> fromNative_(const String& str)
      {
        return fromNative_(str.c_str());
      }

      // Converts from a wide-character string to a narrow-character string.
      inline static String toNative_(const XMLCh* str)
      {
        return String(unique_xerces_ptr<char>(xercesc::XMLString::transcode(str)).get());
      }

      // Converts from a wide-character string to a narrow-character string.
      inline static String toNative_(const unique_xerces_ptr<XMLCh>& str)
      {
        return toNative_(str.get());
      }


public:
      /// Constructor
      StringManager();

      /// Destructor
      ~StringManager();

      /// Transcode the supplied C string to a xerces string
      inline static XercesString convert(const char * str)
      {
        return fromNative_(str).get();
      }

      /// Transcode the supplied C++ string to a xerces string
      inline static XercesString convert(const std::string & str)
      {
        return fromNative_(str.c_str()).get();
      }

      /// Transcode the supplied OpenMS string to a xerces string
      inline static XercesString convert(const String & str)
      {
        return fromNative_(str.c_str()).get();
      }

      /// Transcode the supplied C string to a xerces string pointer
      inline static unique_xerces_ptr<XMLCh> convertPtr(const char * str)
      {
        return fromNative_(str);
      }

      /// Transcode the supplied C++ string to a xerces string pointer
      inline static unique_xerces_ptr<XMLCh> convertPtr(const std::string & str)
      {
        return fromNative_(str.c_str());
      }

      /// Transcode the supplied OpenMS string to a xerces string pointer
      inline static unique_xerces_ptr<XMLCh> convertPtr(const String & str)
      {
        return fromNative_(str.c_str());
      }

      /// Transcode the supplied XMLCh* to a String
      inline static String convert(const XMLCh * str)
      {
        return toNative_(str);
      }

      /**
       * @brief Transcodes the supplied XMLCh* and appends it to the OpenMS String
       *
       * @note Assumes that the XMLCh* only contains ASCII characters
       *
      */
      static void appendASCII(const XMLCh * str, const XMLSize_t length, String & result);

    };

    /**
        @brief Base class for XML handlers.
    */
    class OPENMS_DLLAPI XMLHandler :
      public xercesc::DefaultHandler
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
      XMLHandler(const String & filename, const String & version);
      /// Destructor
      ~XMLHandler() override;

      /// Release internal memory used for parsing (call
      void reset();


      /**
          @name Reimplemented XERCES-C error handlers

          These methods forward the error message to our own error handlers below.
      */
      //@{
      void fatalError(const xercesc::SAXParseException & exception) override;
      void error(const xercesc::SAXParseException & exception) override;
      void warning(const xercesc::SAXParseException & exception) override;
      //@}

      /// Fatal error handler. Throws a ParseError exception
      void fatalError(ActionMode mode, const String & msg, UInt line = 0, UInt column = 0) const;
      /// Error handler for recoverable errors.
      void error(ActionMode mode, const String & msg, UInt line = 0, UInt column = 0) const;
      /// Warning handler.
      void warning(ActionMode mode, const String & msg, UInt line = 0, UInt column = 0) const;

      /// Parsing method for character data
      void characters(const XMLCh * const chars, const XMLSize_t length) override;
      /// Parsing method for opening tags
      void startElement(const XMLCh * const uri, const XMLCh * const localname, const XMLCh * const qname, const xercesc::Attributes & attrs) override;
      /// Parsing method for closing tags
      void endElement(const XMLCh * const uri, const XMLCh * const localname, const XMLCh * const qname) override;

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
      static String writeXMLEscape(const String& to_escape)
      {
        String _copy = to_escape;
        // has() is cheap, so check before calling substitute(), since substitute() will usually happen rarely
        if (_copy.has('&')) _copy.substitute("&","&amp;");
        if (_copy.has('>')) _copy.substitute(">","&gt;");
        if (_copy.has('"')) _copy.substitute("\"","&quot;");
        if (_copy.has('<')) _copy.substitute("<","&lt;");
        if (_copy.has('\'')) _copy.substitute("'","&apos;");

        return _copy;
      }

      /**
      *  @brief Convert an XSD type (e.g. 'xsd:double') to a DataValue.
      * 
      *  Not all conversions are supported yet, due to DataValue using an Int64 as the largest possible integer.
      *  Thus, if the @p value contains a large UInt64, conversion will fail.
      *  Value ranges are currently also not checked, only for XSD types which happen to match the internal representation.
      * 
      *  @param type An XSD type. If the type is not supported, the returned type will be a string
      *  @param value The value in sting format, e.g. "123.34"
      *  @return The Datavalue with the respective type (double, int or string)
      *  @throws Exception::ConversionError if the value does not fit into the internal representation or (for few types) exceeds the XSD specs.
      * 
      */
      static DataValue fromXSDString(const String& type, const String& value)
      {
        DataValue data_value;
        // float type
        if (type == "xsd:double" || type == "xsd:float" || type == "xsd:decimal")
        {
          data_value = DataValue(value.toDouble());
        }
        // <=32 bit integer types
        else if (type == "xsd:byte" ||          // 8bit signed
                 type == "xsd:int" ||           // 32bit signed
                 type == "xsd:unsignedShort" || // 16bit unsigned
                 type == "xsd:short" ||         // 16bit signed
                 type == "xsd:unsignedByte" || type == "xsd:unsignedInt")
        {
          data_value = DataValue(value.toInt32());
        }
        // 64 bit integer types
        else if (type == "xsd:long" || type == "xsd:unsignedLong" ||       // 64bit signed or unsigned respectively
                 type == "xsd:integer" || type == "xsd:negativeInteger" || // any 'integer' has arbitrary size... but we have to cope with 64bit for now.
                 type == "xsd:nonNegativeInteger" || type == "xsd:nonPositiveInteger" || type == "xsd:positiveInteger")
        {
          data_value = DataValue(value.toInt64()); // internally a signed 64-bit integer. So if someone uses 2^64-1 as value, toInt64() will raise an exception...
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

         @param cv A CV, usually the PSI-MS CV, see ControlledVocabulary::getPSIMSCV()
         @param parent_tag The tag which encloses the \<cvParam\>
         @param accession The accession from the 'accession' attribute of the \<cvParam\>
         @param name The name from the 'name' attribute of the \<cvParam\>
         @param value The value from the 'value' attribute of the \<cvParam\>
         @param unit_accession The unit_accession from the 'unitAccession' attribute of the \<cvParam\>
         @return DataValue::EMPTY if a conversion error occured (e.g. if @p value could not be converted to an integer for an @p accession which requires an integer) or the DataValue upon success
      */
      DataValue cvParamToValue(const ControlledVocabulary& cv, const String& parent_tag, 
                               const String& accession, const String& name, const String& value,
                               const String& unit_accession) const;

      /**
         @brief Convert the value of a <em>\<cvParam value=.\></em> (as commonly found in PSI schemata) to the DataValue with the correct type (e.g. int) according to
                the type stored in the CV (usually PSI-MS CV), as well as set its unit.

         @param cv A CV, usually the PSI-MS CV, see ControlledVocabulary::getPSIMSCV()
         @param raw_term Represenation of the raw data (i.e. all strings) from a \<cvParam ...\> without the conversion to a specific value type
         @return DataValue::EMPTY if a conversion error occured (e.g. if @p value could not be converted to an integer for an @p accession which requires an integer) or the DataValue upon success
      */
      DataValue cvParamToValue(const ControlledVocabulary& cv, const CVTerm& raw_term) const;

      /// throws a ParseError if protIDs are not unique, i.e. PeptideIDs will be randomly assigned (bad!)
      /// Should be called before writing any ProtIDs to file
      void checkUniqueIdentifiers_(const std::vector<ProteinIdentification>& prot_ids) const;

protected:
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

      /// parse only until total number of scans and chroms have been determined from attributes
      LOADDETAIL load_detail_; 


      /// Returns if two Xerces strings are equal
      inline bool equal_(const XMLCh * a, const XMLCh * b) const
      {
        return xercesc::XMLString::compareString(a, b) == 0;
      }

      ///@name General MetaInfo handling (for idXML, featureXML, consensusXML)
      //@{

      /// Writes the content of MetaInfoInterface to the file
      void writeUserParam_(const String & tag_name, std::ostream & os, const MetaInfoInterface & meta, UInt indent) const;

      //@}

      ///@name controlled vocabulary handling methods
      //@{

      /// Array of CV term lists (one sublist denotes one term and it's children)
      std::vector<std::vector<String> > cv_terms_;

      /// Converts @p term to the index of the term in the cv_terms_ entry @p section
      /// If the term is not found, @p result_on_error is returned (0 by default)
      SignedSize cvStringToEnum_(const Size section, const String & term, const char * message, const SignedSize result_on_error = 0);

      //@}

      ///@name String conversion
      //@{

      /// Conversion of a String to an integer value
      inline Int asInt_(const String & in) const
      {
        Int res = 0;
        try
        {
          res = in.toInt();
        }
        catch (Exception::ConversionError&)
        {
          error(LOAD, String("Int conversion error of \"") + in + "\"");
        }
        return res;
      }

      /// Conversion of a Xerces string to an integer value
      inline Int asInt_(const XMLCh * in) const
      {
        return xercesc::XMLString::parseInt(in);
      }

      /// Conversion of a String to an unsigned integer value
      inline UInt asUInt_(const String & in) const
      {
        UInt res = 0;
        try
        {
          Int tmp = in.toInt();
          if (tmp < 0)
          {
            throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "");
          }
          res = UInt(tmp);
        }
        catch (Exception::ConversionError& )
        {
          error(LOAD, String("UInt conversion error of \"") + in + "\"");
        }
        return res;
      }

      /// Conversion of a String to a double value
      inline double asDouble_(const String & in) const
      {
        double res = 0.0;
        try
        {
          res = in.toDouble();
        }
        catch (Exception::ConversionError& )
        {
          error(LOAD, String("Double conversion error of \"") + in + "\"");
        }
        return res;
      }

      /// Conversion of a String to a float value
      inline float asFloat_(const String & in) const
      {
        float res = 0.0;
        try
        {
          res = in.toFloat();
        }
        catch (Exception::ConversionError& )
        {
          error(LOAD, String("Float conversion error of \"") + in + "\"");
        }
        return res;
      }

      /**
          @brief Conversion of a string to a boolean value

          'true', 'false', '1' and '0' are accepted.

          @n For all other values a parse error is produced.
      */
      inline bool asBool_(const String & in) const
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
          error(LOAD, String("Boolean conversion error of \"") + in + "\"");
        }
        return false;
      }

      /// Conversion of a xs:datetime string to a DateTime value
      inline DateTime asDateTime_(String date_string) const
      {
        DateTime date_time;
        if (!date_string.empty())
        {
          try
          {
            //strip away milliseconds
            date_string.trim();
            date_string = date_string.substr(0, 19);
            date_time.set(date_string);
          }
          catch (Exception::ParseError& /*err*/ )
          {
            error(LOAD, String("DateTime conversion error of \"") + date_string + "\"");
          }
        }
        return date_time;
      }

      //@}

      ///@name Accessing attributes
      //@{

      /// Converts an attribute to a String
      inline String attributeAsString_(const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convertPtr(name).get());
        if (val == nullptr) fatalError(LOAD, String("Required attribute '") + name + "' not present!");
        return sm_.convert(val);
      }

      /// Converts an attribute to a Int
      inline Int attributeAsInt_(const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convertPtr(name).get());
        if (val == nullptr) fatalError(LOAD, String("Required attribute '") + name + "' not present!");
        return xercesc::XMLString::parseInt(val);
      }

      /// Converts an attribute to a double
      inline double attributeAsDouble_(const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convertPtr(name).get());
        if (val == nullptr) fatalError(LOAD, String("Required attribute '") + name + "' not present!");
        return String(sm_.convert(val)).toDouble();
      }

      /// Converts an attribute to a DoubleList
      inline DoubleList attributeAsDoubleList_(const xercesc::Attributes & a, const char * name) const
      {
        String tmp(expectList_(attributeAsString_(a, name)));
        return ListUtils::create<double>(tmp.substr(1, tmp.size() - 2));
      }

      /// Converts an attribute to an IntList
      inline IntList attributeAsIntList_(const xercesc::Attributes & a, const char * name) const
      {
        String tmp(expectList_(attributeAsString_(a, name)));
        return ListUtils::create<Int>(tmp.substr(1, tmp.size() - 2));
      }

      /// Converts an attribute to an StringList
      inline StringList attributeAsStringList_(const xercesc::Attributes & a, const char * name) const
      {
        String tmp(expectList_(attributeAsString_(a, name)));         
        StringList tmp_list = ListUtils::create<String>(tmp.substr(1, tmp.size() - 2)); // between [ and ]
  
        if (tmp.hasSubstring("\\|")) // check full string for escaped comma
        {
          for (String& s : tmp_list)
          {
            s.substitute("\\|", ",");
          }          
        }
        return tmp_list;
      }

      /**
          @brief Assigns the attribute content to the String @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsString_(String & value, const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convertPtr(name).get());
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
      inline bool optionalAttributeAsInt_(Int & value, const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convertPtr(name).get());
        if (val != nullptr)
        {
          value = xercesc::XMLString::parseInt(val);
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the UInt @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsUInt_(UInt & value, const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convertPtr(name).get());
        if (val != nullptr)
        {
          value = xercesc::XMLString::parseInt(val);
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the double @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsDouble_(double & value, const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convertPtr(name).get());
        if (val != nullptr)
        {
          value = String(sm_.convert(val)).toDouble();
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the DoubleList @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsDoubleList_(DoubleList & value, const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convertPtr(name).get());
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
      inline bool optionalAttributeAsStringList_(StringList & value, const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convertPtr(name).get());
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
      inline bool optionalAttributeAsIntList_(IntList & value, const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convertPtr(name).get());
        if (val != nullptr)
        {
          value = attributeAsIntList_(a, name);
          return true;
        }
        return false;
      }

      /// Converts an attribute to a String
      inline String attributeAsString_(const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val == nullptr) fatalError(LOAD, String("Required attribute '") + sm_.convert(name) + "' not present!");
        return sm_.convert(val);
      }

      /// Converts an attribute to a Int
      inline Int attributeAsInt_(const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val == nullptr) fatalError(LOAD, String("Required attribute '") + sm_.convert(name) + "' not present!");
        return xercesc::XMLString::parseInt(val);
      }

      /// Converts an attribute to a double
      inline double attributeAsDouble_(const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val == nullptr) fatalError(LOAD, String("Required attribute '") + sm_.convert(name) + "' not present!");
        return sm_.convert(val).toDouble();
      }

      /// Converts an attribute to a DoubleList
      inline DoubleList attributeAsDoubleList_(const xercesc::Attributes & a, const XMLCh * name) const
      {
        String tmp(expectList_(attributeAsString_(a, name)));
        return ListUtils::create<double>(tmp.substr(1, tmp.size() - 2));
      }

      /// Converts an attribute to a IntList
      inline IntList attributeAsIntList_(const xercesc::Attributes & a, const XMLCh * name) const
      {
        String tmp(expectList_(attributeAsString_(a, name)));
        return ListUtils::create<Int>(tmp.substr(1, tmp.size() - 2));
      }

      /// Converts an attribute to a StringList
      inline StringList attributeAsStringList_(const xercesc::Attributes & a, const XMLCh * name) const
      {
        String tmp(expectList_(attributeAsString_(a, name)));
        StringList tmp_list = ListUtils::create<String>(tmp.substr(1, tmp.size() - 2)); // between [ and ]

        if (tmp.hasSubstring("\\|")) // check full string for escaped comma
        {
          for (String& s : tmp_list)
          {
            s.substitute("\\|", ",");
          }          
        }
        return tmp_list;
      }

      /// Assigns the attribute content to the String @a value if the attribute is present
      inline bool optionalAttributeAsString_(String& value, const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val != nullptr)
        {
          value = sm_.convert(val);
          return !value.empty();
        }
        return false;
      }

      /// Assigns the attribute content to the Int @a value if the attribute is present
      inline bool optionalAttributeAsInt_(Int & value, const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val != nullptr)
        {
          value = xercesc::XMLString::parseInt(val);
          return true;
        }
        return false;
      }

      /// Assigns the attribute content to the UInt @a value if the attribute is present
      inline bool optionalAttributeAsUInt_(UInt & value, const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val != nullptr)
        {
          value = xercesc::XMLString::parseInt(val);
          return true;
        }
        return false;
      }

      /// Assigns the attribute content to the double @a value if the attribute is present
      inline bool optionalAttributeAsDouble_(double & value, const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val != nullptr)
        {
          value = sm_.convert(val).toDouble();
          return true;
        }
        return false;
      }

      /**
          @brief Assigns the attribute content to the DoubleList @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsDoubleList_(DoubleList & value, const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
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
      inline bool optionalAttributeAsIntList_(IntList & value, const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
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
      inline bool optionalAttributeAsStringList_(StringList & value, const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
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

      inline const String& expectList_(const String& str) const
      {
        if (!(str.hasPrefix('[') && str.hasSuffix(']')))
        {
          fatalError(LOAD, String("List argument is not a string representation of a list!"));
        }
        return str;
      }

    };

  }   // namespace Internal
} // namespace OpenMS


