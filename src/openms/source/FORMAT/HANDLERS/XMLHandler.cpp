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

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <iostream>
#include <vector>
#include <string>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
  namespace Internal
  {

    XMLHandler::XMLHandler(const String & filename, const String & version) :
      error_message_(""),
      file_(filename),
      version_(version)
    {
    }

    XMLHandler::~XMLHandler()
    {
    }

    void XMLHandler::reset()
    {
    }

    void XMLHandler::fatalError(const SAXParseException & exception)
    {
      fatalError(LOAD, sm_.convert(exception.getMessage()), exception.getLineNumber(), exception.getColumnNumber());
    }

    void XMLHandler::error(const SAXParseException & exception)
    {
      error(LOAD, sm_.convert(exception.getMessage()), exception.getLineNumber(), exception.getColumnNumber());
    }

    void XMLHandler::warning(const SAXParseException & exception)
    {
      warning(LOAD, sm_.convert(exception.getMessage()), exception.getLineNumber(), exception.getColumnNumber());
    }

    void XMLHandler::fatalError(ActionMode mode, const String & msg, UInt line, UInt column) const
    {
      if (mode == LOAD)
      {
        error_message_ =  String("While loading '") + file_ + "': " + msg;
      }
      else if (mode == STORE)
      {
        error_message_ =  String("While storing '") + file_ + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message_ += String("( in line ") + line + " column " + column + ")";
      }

      // test if file has the wrong extension and is therefore passed to the wrong parser
      FileTypes::Type ft_name = FileHandler::getTypeByFileName(file_);
      FileTypes::Type ft_content = FileHandler::getTypeByContent(file_);
      if (ft_name != ft_content)
      {
        error_message_ += String("\nProbable cause: The file suffix (") + FileTypes::typeToName(ft_name)
                          + ") does not match the file content (" + FileTypes::typeToName(ft_content) + "). "
                          + "Rename the file to fix this.";
      }

      LOG_FATAL_ERROR << error_message_ << std::endl;
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_, error_message_);
    }

    void XMLHandler::error(ActionMode mode, const String & msg, UInt line, UInt column) const
    {
      if (mode == LOAD)
      {
        error_message_ =  String("Non-fatal error while loading '") + file_ + "': " + msg;
      }
      else if (mode == STORE)
      {
        error_message_ =  String("Non-fatal error while storing '") + file_ + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message_ += String("( in line ") + line + " column " + column + ")";
      }
      LOG_ERROR << error_message_ << std::endl;
    }

    void XMLHandler::warning(ActionMode mode, const String & msg, UInt line, UInt column) const
    {
      if (mode == LOAD)
      {
        error_message_ =  String("While loading '") + file_ + "': " + msg;
      }
      else if (mode == STORE)
      {
        error_message_ =  String("While storing '") + file_ + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message_ += String("( in line ") + line + " column " + column + ")";
      }

// warn only in Debug mode but suppress warnings in release mode (more happy users)
#ifdef OPENMS_ASSERTIONS
      LOG_WARN << error_message_ << std::endl;
#else
      LOG_DEBUG << error_message_ << std::endl;
#endif

    }

    void XMLHandler::characters(const XMLCh * const /*chars*/, const XMLSize_t /*length*/)
    {
    }

    void XMLHandler::startElement(const XMLCh * const /*uri*/, const XMLCh * const /*localname*/, const XMLCh * const /*qname*/, const Attributes & /*attrs*/)
    {
    }

    void XMLHandler::endElement(const XMLCh * const /*uri*/, const XMLCh * const /*localname*/, const XMLCh * const /*qname*/)
    {
    }

    void XMLHandler::writeTo(std::ostream & /*os*/)
    {
    }

    String XMLHandler::errorString()
    {
      return error_message_;
    }

    void XMLHandler::writeUserParam_(const String & tag_name, std::ostream & os, const MetaInfoInterface & meta, UInt indent) const
    {
      std::vector<String> keys;
      meta.getKeys(keys);

      for (Size i = 0; i != keys.size(); ++i)
      {
        os << String(indent, '\t') << "<" << writeXMLEscape(tag_name) << " type=\"";

        DataValue d = meta.getMetaValue(keys[i]);
        String val;
        // determine type
        if (d.valueType() == DataValue::STRING_VALUE || d.valueType() == DataValue::EMPTY_VALUE)
        {
          os << "string";
          val = d;
        }
        else if (d.valueType() == DataValue::INT_VALUE)
        {
          os << "int";
          val = d;
        }
        else if (d.valueType() == DataValue::DOUBLE_VALUE)
        {
          os << "float";
          val = d;
        }
        else if (d.valueType() == DataValue::INT_LIST)
        {
          os << "intList";
          val = d.toString();
        }
        else if (d.valueType() == DataValue::DOUBLE_LIST)
        {
          os << "floatList";
          val = d.toString();
        }
        else if (d.valueType() == DataValue::STRING_LIST)
        {
          os << "stringList";
          StringList sld = d;
          // concatenate manually, as operator<< inserts spaces, which are bad
          // for reconstructing the list
          val = "[" + ListUtils::concatenate(sld, ",") + "]";
        }
        else
        {
          throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
        os << "\" name=\"" << keys[i] << "\" value=\"" << writeXMLEscape(val) << "\"/>" << "\n";
      }
    }

    //*******************************************************************************************************************
    
    StringManager::StringManager()
    {
    }

    StringManager::~StringManager()
    {
    }

    void StringManager::appendASCII(const XMLCh * chars, const XMLSize_t length, String & result)
    {
      // XMLCh are characters in UTF16 (usually stored as 16bit unsigned
      // short but this is not guaranteed).
      // We know that the Base64 string here can only contain plain ASCII
      // and all bytes except the least significant one will be zero. Thus
      // we can convert to char directly (only keeping the least
      // significant byte).
      const XMLCh* it = chars;
      const XMLCh* end = it + length;

      size_t curr_size = result.size();
      result.resize(curr_size + length);
      std::string::iterator str_it = result.begin();
      std::advance(str_it, curr_size);
      while (it!=end)
      {   
        *str_it = (char)*it;
        ++str_it;
        ++it;
      }

      // This is ca. 50 % faster than 
      // for (size_t i = 0; i < length; i++)
      // {
      //   result[curr_size + i] = (char)chars[i];
      // }

    }

  }   // namespace Internal

} // namespace OpenMS
