// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_XMLHANDLER_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h> // StringList
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <xercesc/util/XMLString.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax/Locator.hpp>
#include <xercesc/sax2/Attributes.hpp>

#include <algorithm>
#include <iosfwd>
#include <string>

namespace OpenMS
{
  namespace Internal
  {

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
      inline XercesString fromNative_(const char* str) const
      {
        XMLCh* ptr(xercesc::XMLString::transcode(str));
        XercesString result(ptr);
        xercesc::XMLString::release(&ptr);
        return result;
      }

      // Converts from a narrow-character string to a wide-charactr string.
      inline XercesString fromNative_(const String& str) const
      {
        return fromNative_(str.c_str());
      }

      // Converts from a wide-character string to a narrow-character string.
      inline String toNative_(const XMLCh* str) const
      {
        char* ptr(xercesc::XMLString::transcode(str));
        String result(ptr);
        xercesc::XMLString::release(&ptr);
        return result;
      }

      // Converts from a wide-character string to a narrow-character string.
      inline String toNative_(const XercesString& str) const
      {
        return toNative_(str.c_str());
      }


public:
      /// Constructor
      StringManager();

      /// Destructor
      ~StringManager();

      /// Transcode the supplied C string to a xerces string
      inline XercesString convert(const char * str) const
      {
        return fromNative_(str);
      }

      /// Transcode the supplied C++ string to a xerces string
      inline XercesString convert(const std::string & str) const
      {
        return fromNative_(str.c_str());
      }

      /// Transcode the supplied OpenMS string to a xerces string
      inline XercesString convert(const String & str) const
      {
        return fromNative_(str.c_str());
      }

      /// Transcode the supplied XMLCh* to a String
      inline String convert(const XMLCh * str) const
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

      /// Default constructor
      XMLHandler(const String & filename, const String & version);
      /// Destructor
      virtual ~XMLHandler();

      /// Release internal memory used for parsing (call
      void reset();


      /**
          @name Reimplemented XERCES-C error handlers

          These methods forward the error message to our own error handlers below.
      */
      //@{
      void fatalError(const xercesc::SAXParseException & exception);
      void error(const xercesc::SAXParseException & exception);
      void warning(const xercesc::SAXParseException & exception);
      //@}

      /// Fatal error handler. Throws a ParseError exception
      void fatalError(ActionMode mode, const String & msg, UInt line = 0, UInt column = 0) const;
      /// Error handler for recoverable errors.
      void error(ActionMode mode, const String & msg, UInt line = 0, UInt column = 0) const;
      /// Warning handler.
      void warning(ActionMode mode, const String & msg, UInt line = 0, UInt column = 0) const;

      /// Parsing method for character data
      virtual void characters(const XMLCh * const chars, const XMLSize_t length);
      /// Parsing method for opening tags
      virtual void startElement(const XMLCh * const uri, const XMLCh * const localname, const XMLCh * const qname, const xercesc::Attributes & attrs);
      /// Parsing method for closing tags
      virtual void endElement(const XMLCh * const uri, const XMLCh * const localname, const XMLCh * const qname);

      /// Writes the contents to a stream.
      virtual void writeTo(std::ostream & /*os*/);

      /// Returns the last error description
      String errorString();

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
      inline SignedSize cvStringToEnum_(const Size section, const String & term, const char * message, const SignedSize result_on_error = 0)
      {
        OPENMS_PRECONDITION(section < cv_terms_.size(), "cvStringToEnum_: Index overflow (section number too large)");

        std::vector<String>::const_iterator it = std::find(cv_terms_[section].begin(), cv_terms_[section].end(), term);
        if (it != cv_terms_[section].end())
        {
          return it - cv_terms_[section].begin();
        }
        else
        {
          warning(LOAD, String("Unexpected CV entry '") + message + "'='" + term + "'");
          return result_on_error;
        }
      }

      //@}

      ///@name String conversion
      //@{

      /// Conversion of a String to an integer value
      inline Int asInt_(const String & in)
      {
        Int res = 0;
        try
        {
          res = in.toInt();
        }
        catch (Exception::ConversionError)
        {
          error(LOAD, String("Int conversion error of \"") + in + "\"");
        }
        return res;
      }

      /// Conversion of a Xerces string to an integer value
      inline Int asInt_(const XMLCh * in)
      {
        return xercesc::XMLString::parseInt(in);
      }

      /// Conversion of a String to an unsigned integer value
      inline UInt asUInt_(const String & in)
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
      inline double asDouble_(const String & in)
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
      inline float asFloat_(const String & in)
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
      inline bool asBool_(const String & in)
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
      inline DateTime asDateTime_(String date_string)
      {
        DateTime date_time;
        if (date_string != "")
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
        const XMLCh * val = a.getValue(sm_.convert(name).c_str());
        if (val == 0) fatalError(LOAD, String("Required attribute '") + name + "' not present!");
        return sm_.convert(val);
      }

      /// Converts an attribute to a Int
      inline Int attributeAsInt_(const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convert(name).c_str());
        if (val == 0) fatalError(LOAD, String("Required attribute '") + name + "' not present!");
        return xercesc::XMLString::parseInt(val);
      }

      /// Converts an attribute to a double
      inline double attributeAsDouble_(const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convert(name).c_str());
        if (val == 0) fatalError(LOAD, String("Required attribute '") + name + "' not present!");
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
        return ListUtils::create<String>(tmp.substr(1, tmp.size() - 2));
      }

      /**
          @brief Assigns the attribute content to the String @a value if the attribute is present

          @return if the attribute was present
      */
      inline bool optionalAttributeAsString_(String & value, const xercesc::Attributes & a, const char * name) const
      {
        const XMLCh * val = a.getValue(sm_.convert(name).c_str());
        if (val != 0)
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
        const XMLCh * val = a.getValue(sm_.convert(name).c_str());
        if (val != 0)
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
        const XMLCh * val = a.getValue(sm_.convert(name).c_str());
        if (val != 0)
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
        const XMLCh * val = a.getValue(sm_.convert(name).c_str());
        if (val != 0)
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
        const XMLCh * val = a.getValue(sm_.convert(name).c_str());
        if (val != 0)
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
        const XMLCh * val = a.getValue(sm_.convert(name).c_str());
        if (val != 0)
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
        const XMLCh * val = a.getValue(sm_.convert(name).c_str());
        if (val != 0)
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
        if (val == 0) fatalError(LOAD, String("Required attribute '") + sm_.convert(name) + "' not present!");
        return sm_.convert(val);
      }

      /// Converts an attribute to a Int
      inline Int attributeAsInt_(const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val == 0) fatalError(LOAD, String("Required attribute '") + sm_.convert(name) + "' not present!");
        return xercesc::XMLString::parseInt(val);
      }

      /// Converts an attribute to a double
      inline double attributeAsDouble_(const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val == 0) fatalError(LOAD, String("Required attribute '") + sm_.convert(name) + "' not present!");
        return String(sm_.convert(val)).toDouble();
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
        return ListUtils::create<String>(tmp.substr(1, tmp.size() - 2));
      }

      /// Assigns the attribute content to the String @a value if the attribute is present
      inline bool optionalAttributeAsString_(String & value, const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val != 0)
        {
          String tmp2 = sm_.convert(val);
          if (tmp2 != "")
          {
            value = tmp2;
            return true;
          }
        }
        return false;
      }

      /// Assigns the attribute content to the Int @a value if the attribute is present
      inline bool optionalAttributeAsInt_(Int & value, const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val != 0)
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
        if (val != 0)
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
        if (val != 0)
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
      inline bool optionalAttributeAsDoubleList_(DoubleList & value, const xercesc::Attributes & a, const XMLCh * name) const
      {
        const XMLCh * val = a.getValue(name);
        if (val != 0)
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
        if (val != 0)
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
        if (val != 0)
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

      inline String expectList_(const String& str) const
      {
        String tmp(str);
        if (!(tmp.hasPrefix('[') && tmp.hasSuffix(']')))
        {
          fatalError(LOAD, String("List argument is not a string representation of a list!"));
        }
        return tmp;
      }

    };

  }   // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_XMLHANDLER_H

