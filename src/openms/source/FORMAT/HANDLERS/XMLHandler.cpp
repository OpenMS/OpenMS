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
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/SIMDe.h>

#include <algorithm>
#include <bit>
#include <set>

using namespace std;

namespace OpenMS::Internal
{

    XMLHandler::XMLHandler(const std::string & filename, const std::string & version) :
      file_(filename),
      version_(version),
      load_detail_(LD_ALLDATA)
    {
    }

    XMLHandler::~XMLHandler() = default;

    void XMLHandler::reset()
    {
    }

    void XMLHandler::fatalError(ActionMode mode, const std::string & msg, UInt line, UInt column) const
    {
      std::string error_message;
      if (mode == LOAD)
      {
        error_message =std::string("While loading '") + file_ + "': " + msg;
	      // test if file has the wrong extension and is therefore passed to the wrong parser
        // only makes sense if we are loading/parsing a file
	      FileTypes::Type ft_name = FileHandler::getTypeByFileName(file_);
        FileTypes::Type ft_content = FileHandler::getTypeByContent(file_);
        if (ft_name != ft_content)
        {
          error_message +=std::string("\nProbable cause: The file suffix (") + FileTypes::typeToName(ft_name)
                          + ") does not match the file content (" + FileTypes::typeToName(ft_content) + "). "
                          + "Rename the file to fix this.";
        }
      }
      else if (mode == STORE)
      {
        error_message =std::string("While storing '") + file_ + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message +=std::string("( in line ") + line + " column " + column + ")";
      }

      OPENMS_LOG_FATAL_ERROR << error_message << std::endl;
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_, error_message);
    }

    void XMLHandler::error(ActionMode mode, const std::string & msg, UInt line, UInt column) const
    {
      std::string error_message;
      if (mode == LOAD)
      {
        error_message =std::string("Non-fatal error while loading '") + file_ + "': " + msg;
      }
      else if (mode == STORE)
      {
        error_message =std::string("Non-fatal error while storing '") + file_ + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message +=std::string("( in line ") + line + " column " + column + ")";
      }
      OPENMS_LOG_ERROR << error_message << std::endl;
    }

    void XMLHandler::warning(ActionMode mode, const std::string & msg, UInt line, UInt column) const
    {
      std::string error_message;
      if (mode == LOAD)
      {
        error_message =std::string("While loading '") + file_ + "': " + msg;
      }
      else if (mode == STORE)
      {
        error_message =std::string("While storing '") + file_ + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message +=std::string("( in line ") + line + " column " + column + ")";
      }

// warn only in Debug mode but suppress warnings in release mode (more happy users)
#ifdef OPENMS_ASSERTIONS
      OPENMS_LOG_WARN << error_message << std::endl;
#else
      OPENMS_LOG_DEBUG << error_message << std::endl;
#endif

    }

    // Native default implementations (overridden by migrated subclasses).
    void XMLHandler::onStartElement(const char16_t * /*qname*/, const XMLAttributes & /*attributes*/)
    {
    }

    void XMLHandler::onEndElement(const char16_t * /*qname*/)
    {
    }

    void XMLHandler::onCharacters(const char16_t * /*chars*/, Size /*length*/)
    {
    }

    void XMLHandler::writeTo(std::ostream & /*os*/)
    {
    }

    SignedSize XMLHandler::cvStringToEnum_(const Size section, const std::string & term, const char * message, const SignedSize result_on_error)
    {
      OPENMS_PRECONDITION(section < cv_terms_.size(), "cvStringToEnum_: Index overflow (section number too large)");

      std::vector<std::string>::const_iterator it = std::find(cv_terms_[section].begin(), cv_terms_[section].end(), term);
      if (it != cv_terms_[section].end())
      {
        return it - cv_terms_[section].begin();
      }
      else
      {
        warning(LOAD,std::string("Unexpected CV entry '") + message + "'='" + term + "'");
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


    DataValue XMLHandler::cvParamToValue(const ControlledVocabulary& cv, const std::string& parent_tag, const std::string& accession, const std::string& name,
                                         const std::string& value, const std::string& unit_accession) const
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
          warning(LOAD,std::string("Name of CV term not correct: '") + term.id + " - " + name + "' should be '" + term.name + "'");
        }
        if (term.obsolete)
        {
          warning(LOAD,std::string("Obsolete CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "'.");
        }
        // values used in wrong places and wrong value types
        if (!value.empty())
        {
          if (term.xref_type == ControlledVocabulary::CVTerm::XRefType::NONE)
          {
            // Quality CV does not state value type :(
            if (!StringUtils::hasPrefix(accession, "PATO:"))
            {
              warning(LOAD,std::string("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' must not have a value. The value is '" + value + "'.");
            }
          }
          else
          {
            switch (term.xref_type)
            {
              // string value can be anything
              case ControlledVocabulary::CVTerm::XRefType::XSD_STRING:
                break;

              // int value => try casting
              case ControlledVocabulary::CVTerm::XRefType::XSD_INTEGER:
              case ControlledVocabulary::CVTerm::XRefType::XSD_NEGATIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XRefType::XSD_POSITIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XRefType::XSD_NON_NEGATIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XRefType::XSD_NON_POSITIVE_INTEGER:
                try
                {
                  cv_value = StringUtils::toInt32(value);
                }
                catch (Exception::ConversionError&)
                {
                  warning(LOAD,std::string("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' must have an integer value. The value is '" + value + "'.");
                  return DataValue::EMPTY;
                }
                break;

              // double value => try casting
              case ControlledVocabulary::CVTerm::XRefType::XSD_DECIMAL:
                try
                {
                  cv_value = StringUtils::toDouble(value);
                }
                catch (Exception::ConversionError&)
                {
                  warning(LOAD,std::string("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' must have a floating-point value. The value is '" + value + "'.");
                  return DataValue::EMPTY;
                }
                break;

              // date string => try conversion
              case ControlledVocabulary::CVTerm::XRefType::XSD_DATE:
                try
                {
                  DateTime tmp;
                  tmp.set(value);
                }
                catch (Exception::ParseError&)
                {
                  warning(LOAD,std::string("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' must be a valid date. The value is '" + value + "'.");
                  return DataValue::EMPTY;
                }
                break;
              
              case ControlledVocabulary::CVTerm::XRefType::XSD_BOOLEAN:
                try
                {
                  cv_value = StringUtils::toLowered(value);
                  cv_value.toBool(); // only works if 'true' or 'false'
                }
                catch (Exception::ConversionError&)
                {
                  warning(LOAD,std::string("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' must have a boolean value (true/false). The value is '" + value + "'.");
                  return DataValue::EMPTY;
                }
                break;

              default:
                warning(LOAD,std::string("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' has the unknown value type '" +
                                ControlledVocabulary::CVTerm::getXRefTypeName(term.xref_type) + "'.");
                break;
            }
          }
        }
        // no value, although there should be a numerical value
        else if (term.xref_type != ControlledVocabulary::CVTerm::XRefType::NONE && term.xref_type != ControlledVocabulary::CVTerm::XRefType::XSD_STRING && // should be numerical
                 !cv.isChildOf(accession, "MS:1000513") // here the value type relates to the binary data array, not the 'value=' attribute!
        )
        {
          warning(LOAD,std::string("The CV term '") + accession + " - " + term.name + "' used in tag '" + parent_tag + "' should have a numerical value. The value is '" + value + "'.");
          return DataValue::EMPTY;
        }
      }
      catch (const Exception::InvalidValue& /*e*/)
      {
        // in 'sample' several external CVs are used (Brenda, GO, ...). Do not warn then.
        if (parent_tag != "sample")
        {
          warning(LOAD,std::string("Unknown cvParam '") + accession + "' in tag '" + parent_tag + "'.");
          return DataValue::EMPTY;
        }
      }

      if (!unit_accession.empty())
      {
        if (StringUtils::hasPrefix(unit_accession, "UO:"))
        {
          cv_value.setUnit(StringUtils::toInt32(StringUtils::suffix(unit_accession, unit_accession.size() - 3)));
          cv_value.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
        }
        else if (StringUtils::hasPrefix(unit_accession, "MS:"))
        {
          cv_value.setUnit(StringUtils::toInt32(StringUtils::suffix(unit_accession, unit_accession.size() - 3)));
          cv_value.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
        }
        else
        {
          warning(LOAD,std::string("Unhandled unit '") + unit_accession + "' in tag '" + parent_tag + "'.");
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
      std::set<std::string> s;
      for (const auto& p : prot_ids)
      {
        if (s.insert(p.getIdentifier()).second == false) // element already existed
        {
          fatalError(ActionMode::STORE, "ProteinIdentification run identifiers are not unique. This can lead to loss of unique PeptideIdentification assignment. Duplicated Protein-ID is:" +
                                        p.getIdentifier());
        }
      }
    }

    void XMLHandler::writeUserParam_(const std::string& tag_name, std::ostream& os, const MetaInfoInterface& meta, UInt indent) const
    {
      std::vector<std::string> keys;
      meta.getKeys(keys);

      std::string val;
      std::string p_prefix = std::string(indent, '\t') + "<" + writeXMLEscape(tag_name) + " type=\"";
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
          val = StringUtils::toStr(d);
        }
        else if (d.valueType() == DataValue::DOUBLE_VALUE)
        {
          os << "float";
          val = StringUtils::toStr(d);
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
          for (std::string& s : sl)
          {
            if (StringUtils::has(s, ',')) StringUtils::substitute(s, ",", "\\|");
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


    //*******************************************************************************************************************

    DoubleList XMLHandler::attributeAsDoubleList_(const XMLAttributes & a, const char * name) const
    {
      std::string tmp(expectList_(attributeAsString_(a, name)));
      return ListUtils::create<double>(StringUtils::substr(tmp, 1, tmp.size() - 2));
    }

    IntList XMLHandler::attributeAsIntList_(const XMLAttributes & a, const char * name) const
    {
      std::string tmp(expectList_(attributeAsString_(a, name)));
      return ListUtils::create<Int>(StringUtils::substr(tmp, 1, tmp.size() - 2));
    }

    StringList XMLHandler::attributeAsStringList_(const XMLAttributes & a, const char * name) const
    {
      std::string tmp(expectList_(attributeAsString_(a, name)));
      StringList tmp_list = ListUtils::create<std::string>(StringUtils::substr(tmp, 1, tmp.size() - 2)); // between [ and ]

      if (StringUtils::hasSubstring(tmp, "\\|")) // check full string for escaped comma
      {
        for (std::string& s : tmp_list)
        {
          StringUtils::substitute(s, "\\|", ",");
        }
      }
      return tmp_list;
    }

    DoubleList XMLHandler::attributeAsDoubleList_(const XMLAttributes & a, const char16_t * name) const
    {
      std::string tmp(expectList_(attributeAsString_(a, name)));
      return ListUtils::create<double>(StringUtils::substr(tmp, 1, tmp.size() - 2));
    }

    IntList XMLHandler::attributeAsIntList_(const XMLAttributes & a, const char16_t * name) const
    {
      std::string tmp(expectList_(attributeAsString_(a, name)));
      return ListUtils::create<Int>(StringUtils::substr(tmp, 1, tmp.size() - 2));
    }

    StringList XMLHandler::attributeAsStringList_(const XMLAttributes & a, const char16_t * name) const
    {
      std::string tmp(expectList_(attributeAsString_(a, name)));
      StringList tmp_list = ListUtils::create<std::string>(StringUtils::substr(tmp, 1, tmp.size() - 2)); // between [ and ]

      if (StringUtils::hasSubstring(tmp, "\\|")) // check full string for escaped comma
      {
        for (std::string& s : tmp_list)
        {
          StringUtils::substitute(s, "\\|", ",");
        }
      }
      return tmp_list;
    }


} // namespace OpenMS::Internal
