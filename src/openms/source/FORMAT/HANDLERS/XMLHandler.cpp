// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/SIMDe.h>

#include <algorithm>
#include <set>

using namespace std;
using namespace xercesc;

namespace OpenMS::Internal
{

    // Specializations for character types, released by XMLString::release
    template<> void shared_xerces_ptr<char>::doRelease_(char* item)
    {
      xercesc::XMLString::release(&item);
    }

    template<> void shared_xerces_ptr<XMLCh>::doRelease_(XMLCh* item)
    {
      xercesc::XMLString::release(&item);
    }

    // Specializations for character types, which needs to be
    // released by XMLString::release
    template <>
    void unique_xerces_ptr<char>::doRelease_(char*& item)
    {
      xercesc::XMLString::release(&item);
    }

    template <>
    void unique_xerces_ptr<XMLCh>::doRelease_(XMLCh*& item)
    {
      xercesc::XMLString::release(&item);
    }

    XMLHandler::XMLHandler(const String & filename, const String & version) :
      file_(filename),
      version_(version),
      load_detail_(LD_ALLDATA)
    {
    }

    XMLHandler::~XMLHandler() = default;

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
      String error_message;
      if (mode == LOAD)
      {
        error_message =  String("While loading '") + file_ + "': " + msg;
	      // test if file has the wrong extension and is therefore passed to the wrong parser
        // only makes sense if we are loading/parsing a file
	      FileTypes::Type ft_name = FileHandler::getTypeByFileName(file_);
        FileTypes::Type ft_content = FileHandler::getTypeByContent(file_);
        if (ft_name != ft_content)
        {
          error_message += String("\nProbable cause: The file suffix (") + FileTypes::typeToName(ft_name)
                          + ") does not match the file content (" + FileTypes::typeToName(ft_content) + "). "
                          + "Rename the file to fix this.";
        }
      }
      else if (mode == STORE)
      {
        error_message =  String("While storing '") + file_ + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message += String("( in line ") + line + " column " + column + ")";
      }

      OPENMS_LOG_FATAL_ERROR << error_message << std::endl;
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_, error_message);
    }

    void XMLHandler::error(ActionMode mode, const String & msg, UInt line, UInt column) const
    {
      String error_message;
      if (mode == LOAD)
      {
        error_message =  String("Non-fatal error while loading '") + file_ + "': " + msg;
      }
      else if (mode == STORE)
      {
        error_message =  String("Non-fatal error while storing '") + file_ + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message += String("( in line ") + line + " column " + column + ")";
      }
      OPENMS_LOG_ERROR << error_message << std::endl;
    }

    void XMLHandler::warning(ActionMode mode, const String & msg, UInt line, UInt column) const
    {
      String error_message;
      if (mode == LOAD)
      {
        error_message =  String("While loading '") + file_ + "': " + msg;
      }
      else if (mode == STORE)
      {
        error_message =  String("While storing '") + file_ + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message += String("( in line ") + line + " column " + column + ")";
      }

// warn only in Debug mode but suppress warnings in release mode (more happy users)
#ifdef OPENMS_ASSERTIONS
      OPENMS_LOG_WARN << error_message << std::endl;
#else
      OPENMS_LOG_DEBUG << error_message << std::endl;
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

    SignedSize XMLHandler::cvStringToEnum_(const Size section, const String & term, const char * message, const SignedSize result_on_error)
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

    /// handlers which support partial loading, implement this method
    /// @throws Exception::NotImplemented
    XMLHandler::LOADDETAIL XMLHandler::getLoadDetail() const
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    /// handlers which support partial loading, implement this method
    /// @throws Exception::NotImplemented
    void XMLHandler::setLoadDetail(const LOADDETAIL /*d*/)
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }


    DataValue XMLHandler::cvParamToValue(const ControlledVocabulary& cv, const String& parent_tag, const String& accession, const String& name,
                                         const String& value, const String& unit_accession) const
    {
      // the actual value stored in the CVParam
      // we assume for now that it is a string value, we update the type later on
      DataValue cv_value = value;

      // Abort on unknown terms
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(accession); // throws Exception::InvalidValue if missing

        // check if term name and parsed name match
        if (name != term.name) 
        {
          warning(LOAD, String("Name of CV term not correct: '") + term.id + " - " + name + "' should be '" + term.name + "'");
        }
        if (term.obsolete)
        {
          warning(LOAD, String("Obsolete CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "'.");
        }
        // values used in wrong places and wrong value types
        if (!value.empty())
        {
          if (term.xref_type == ControlledVocabulary::CVTerm::NONE)
          {
            // Quality CV does not state value type :(
            if (!accession.hasPrefix("PATO:"))
            {
              warning(LOAD, String("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' must not have a value. The value is '" + value + "'.");
            }
          }
          else
          {
            switch (term.xref_type)
            {
              // string value can be anything
              case ControlledVocabulary::CVTerm::XSD_STRING:
                break;

              // int value => try casting
              case ControlledVocabulary::CVTerm::XSD_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NEGATIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_POSITIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NON_NEGATIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NON_POSITIVE_INTEGER:
                try
                {
                  cv_value = value.toInt();
                }
                catch (Exception::ConversionError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' must have an integer value. The value is '" + value + "'.");
                  return DataValue::EMPTY;
                }
                break;

              // double value => try casting
              case ControlledVocabulary::CVTerm::XSD_DECIMAL:
                try
                {
                  cv_value = value.toDouble();
                }
                catch (Exception::ConversionError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' must have a floating-point value. The value is '" + value + "'.");
                  return DataValue::EMPTY;
                }
                break;

              // date string => try conversion
              case ControlledVocabulary::CVTerm::XSD_DATE:
                try
                {
                  DateTime tmp;
                  tmp.set(value);
                }
                catch (Exception::ParseError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' must be a valid date. The value is '" + value + "'.");
                  return DataValue::EMPTY;
                }
                break;
              
              case ControlledVocabulary::CVTerm::XSD_BOOLEAN:
                try
                {
                  cv_value = String(value).toLower();
                  cv_value.toBool(); // only works if 'true' or 'false'
                }
                catch (Exception::ConversionError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' must have a boolean value (true/false). The value is '" + value + "'.");
                  return DataValue::EMPTY;
                }
                break;

              default:
                warning(LOAD, String("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' has the unknown value type '" +
                                ControlledVocabulary::CVTerm::getXRefTypeName(term.xref_type) + "'.");
                break;
            }
          }
        }
        // no value, although there should be a numerical value
        else if (term.xref_type != ControlledVocabulary::CVTerm::NONE && term.xref_type != ControlledVocabulary::CVTerm::XSD_STRING && // should be numerical
                 !cv.isChildOf(accession, "MS:1000513") // here the value type relates to the bits data array, not the 'value=' attribute!
        )
        {
          warning(LOAD, String("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' should have a numerical value. The value is '" + value + "'.");
          return DataValue::EMPTY;
        }
      }
      catch (const Exception::InvalidValue& /*e*/)
      {
        // in 'sample' several external CVs are used (Brenda, GO, ...). Do not warn then.
        if (parent_tag != "sample")
        {
          warning(LOAD, String("Unknown cvParam '") + accession + "' in tag '" + parent_tag + "'.");
          return DataValue::EMPTY;
        }
      }

      if (!unit_accession.empty())
      {
        if (unit_accession.hasPrefix("UO:"))
        {
          cv_value.setUnit(unit_accession.suffix(unit_accession.size() - 3).toInt());
          cv_value.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
        }
        else if (unit_accession.hasPrefix("MS:"))
        {
          cv_value.setUnit(unit_accession.suffix(unit_accession.size() - 3).toInt());
          cv_value.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
        }
        else
        {
          warning(LOAD, String("Unhandled unit '") + unit_accession + "' in tag '" + parent_tag + "'.");
        }
      }

      return cv_value;
    }

    DataValue XMLHandler::cvParamToValue(const ControlledVocabulary& cv, const CVTerm& raw_term) const
    {
      return cvParamToValue(cv, "?", raw_term.getAccession(), raw_term.getName(), raw_term.getValue(), raw_term.getUnit().accession);
    }

    void XMLHandler::checkUniqueIdentifiers_(const std::vector<ProteinIdentification>& prot_ids) const
    {
      std::set<String> s;
      for (const auto& p : prot_ids)
      {
        if (s.insert(p.getIdentifier()).second == false) // element already existed
        {
          fatalError(ActionMode::STORE, "ProteinIdentification run identifiers are not unique. This can lead to loss of unique PeptideIdentification assignment. Duplicated Protein-ID is:" +
                                        p.getIdentifier());
        }
      }
    }

    void XMLHandler::writeUserParam_(const String& tag_name, std::ostream& os, const MetaInfoInterface& meta, UInt indent) const
    {
      std::vector<String> keys;
      meta.getKeys(keys);

      String val;
      String p_prefix = String(indent, '\t') + "<" + writeXMLEscape(tag_name) + " type=\"";
      for (Size i = 0; i != keys.size(); ++i)
      {
        os << p_prefix;

        const DataValue& d = meta.getMetaValue(keys[i]);
        // determine type
        if (d.valueType() == DataValue::STRING_VALUE || d.valueType() == DataValue::EMPTY_VALUE)
        {
          os << "string";
          val = writeXMLEscape(d);
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
          // List elements are separated by comma. In the rare case of comma inside individual strings
          // we replace them by an escape symbol '\\|'. 
          // This allows distinguishing commas as element separator and normal string character and reconstruct the list.
          StringList sl = d.toStringList();
          for (String& s : sl)
          {
            if (s.has(',')) s.substitute(",", "\\|");
          }
          val = "[" + writeXMLEscape(ListUtils::concatenate(sl, ",")) + "]";
        }
        else
        {
          throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
        os << "\" name=\"" << keys[i] << "\" value=\"" << val << "\"/>\n";
      }
    }

    size_t StringManager::strLength(const XMLCh* input_ptr) {
      if (input_ptr == nullptr) {
          return 0;
      }
  
      XMLSize_t processed_chars = 0;
      const XMLCh* pos_ptr = input_ptr;
  
      // Verarbeite einzelne Zeichen, bis der Pointer 16-Byte-aligned ist
      uintptr_t ptr_value = reinterpret_cast<uintptr_t>(pos_ptr);
      size_t misalignment = ptr_value & 0xF; // Berechnet Misalignment als (Adresswert) mod 16
      size_t chars_to_align = misalignment ? (16 - misalignment) / sizeof(XMLCh) : 0;
  
      // Vorverarbeitung einzelner Zeichen bis zum Alignment oder bis zum Ende des Strings
      for (size_t i = 0; i < chars_to_align; ++i) {
          if (*pos_ptr == 0) {
              return processed_chars;
          }
          ++pos_ptr;
          ++processed_chars;
      }
  
      // Hauptschleife mit SIMD-Operationen
      while (true) {
          // SIMD-Operation
          simde__m128i bits = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(pos_ptr));
          simde__m128i zero = simde_mm_setzero_si128();
          simde__m128i cmp_zero = simde_mm_cmpeq_epi16(bits, zero);
          uint16_t zero_mask = simde_mm_movemask_epi8(cmp_zero);
  
          if (zero_mask != 0x0000) {
              int byte_pos_zero = __builtin_ctz(zero_mask);
              int char_pos_zero = byte_pos_zero / 2;
              return processed_chars + static_cast<>(char_pos_zero);
          }
  
          // 8 Zeichen (16 Bytes) wurden verarbeitet, keine Null gefunden
          pos_ptr += 8;
          processed_chars += 8;
      }
  
      // Diese Zeile wird nie erreicht
      return processed_chars;
  }

    void StringManager::compress64(const XMLCh* inputIt, char* outputIt)
    {
      simde__m128i bits = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(inputIt));
    
      // Select every second byte (little-endian lower byte of each UTF-16 character)
      const simde__m128i shuffleMask = simde_mm_setr_epi8(
        0, 2, 4, 6, 8, 10, 12, 14,
        -1, -1, -1, -1, -1, -1, -1, -1
      );
    
      simde__m128i compressed = simde_mm_shuffle_epi8(bits, shuffleMask);
    
      // Store the lower 64 bits (8 ASCII characters)
      simde_mm_storel_epi64(reinterpret_cast<simde__m128i*>(outputIt), compressed);
    }

    bool StringManager::isASCII(const XMLCh* chars, const XMLSize_t length)
  {
    if (length == 0)
    {
      return false;
    }

    Size quotient = length / 8;
    Size remainder = length % 8;

    const XMLCh* inputPtr = chars;
    simde__m128i mask = simde_mm_set1_epi16(0xFF00);
    bool bitmask = true;

    // Process blocks of 8 UTF-16 characters using SIMD
    for (Size i = 0; i < quotient; ++i)
    {
      simde__m128i bits = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(inputPtr));
      simde__m128i zero = simde_mm_setzero_si128();
      simde__m128i andOp = simde_mm_and_si128(bits, mask);
      simde__m128i cmp = simde_mm_cmpeq_epi16(andOp, zero);

      if (simde_mm_movemask_epi8(cmp) != 0xFFFF)
      {
        bitmask = false;
        break;
      }

      inputPtr += 8;
    }

  // Check remaining characters individually
    for (Size i = 0; i < remainder && bitmask; ++i)
    {
      if (inputPtr[i] & 0xFF00)
      {
        bitmask = false;
        break;
      }
    }

    return bitmask;
  }

    void StringManager::appendASCII(const XMLCh* chars, const XMLSize_t length, String& result)
    {
      // XMLCh are characters in UTF16 (usually stored as 16-bit unsigned
      // short but this is not guaranteed).
      // We know that the Base64 string here can only contain plain ASCII
      // and all bytes except the least significant one will be zero. Thus
      // we can convert to char directly (only keeping the least
      // significant byte).
    
      Size quotient = length / 8;
      Size remainder = length % 8;
    
      const XMLCh* inputPtr = chars;
    
      Size currentSize = result.size();
      result.resize(currentSize + length);
      char* outputPtr = &result[currentSize];
    
      // Copy blocks of 8 characters at a time
      for (Size i = 0; i < quotient; ++i)
      {
        compress64(inputPtr, outputPtr);
        inputPtr += 8;
        outputPtr += 8;
      }
    
      // Copy any remaining characters individually
      for (Size i = 0; i < remainder; ++i)
      {
        outputPtr[i] = static_cast<char>(inputPtr[i] & 0xFF);
      }
    }

    //*******************************************************************************************************************
    
    StringManager::StringManager()
    = default;

    StringManager::~StringManager()
    = default;

    


} // namespace OpenMS   // namespace Internal
