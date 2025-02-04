// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TransformationXMLFile.h>

#include <fstream>

using namespace std;

namespace OpenMS
{

  TransformationXMLFile::TransformationXMLFile() :
    XMLHandler("", "1.1"),
    XMLFile("/SCHEMAS/TrafoXML_1_1.xsd", "1.1"),
    params_(), data_(), model_type_()
  {
  }

  void TransformationXMLFile::load(const String& filename, TransformationDescription& transformation, bool fit_model)
  {
    //Filename for error messages in XMLHandler
    file_ = filename;

    params_.clear();
    data_.clear();
    model_type_.clear();

    parse_(filename, this);

    transformation.setDataPoints(data_);

    if (fit_model)
    {
      transformation.fitModel(model_type_, params_);
    }
  }

  void TransformationXMLFile::store(const String& filename, const TransformationDescription& transformation)
  {
    if (transformation.getModelType().empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "will not write a transformation with empty name");
    }

    //open stream
    std::ofstream os(filename.c_str());
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    os.precision(writtenDigits<double>(0.0));

    //write header
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    //add XSLT file if it can be found
    os << "<TrafoXML version=\"" << getVersion() << "\" xsi:noNamespaceSchemaLocation=\"https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/"
       << schema_location_.suffix('/') << "\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

    // open tag
    os << "\t<Transformation name=\"" << transformation.getModelType()
       << "\">\n";

    // write parameters
    const Param& params = transformation.getModelParameters();
    for (Param::ParamIterator it = params.begin(); it != params.end(); ++it)
    {
      if (it->value.valueType() != ParamValue::EMPTY_VALUE)
      {
        switch (it->value.valueType())
        {
        case ParamValue::INT_VALUE:
          os << "\t\t<Param  type=\"int\" name=\"" << it->name << "\" value=\"" << it->value.toString() << "\"/>\n";
          break;

        case ParamValue::DOUBLE_VALUE:
          os << "\t\t<Param  type=\"float\" name=\"" << it->name << "\" value=\"" << it->value.toString() << "\"/>\n";
          break;

        case ParamValue::STRING_VALUE:
        case ParamValue::STRING_LIST:
        case ParamValue::INT_LIST:
        case ParamValue::DOUBLE_LIST:
          os << "\t\t<Param  type=\"string\" name=\"" << it->name << "\" value=\"" << it->value.toString() << "\"/>\n";
          break;

        default:         // no other value types are supported!
          fatalError(STORE, String("Unsupported parameter type of parameter '") + it->name + "'");
          break;
        }
      }
    }

    //write pairs
    if (!transformation.getDataPoints().empty())
    {
      os << "\t\t<Pairs count=\"" << transformation.getDataPoints().size() << "\">\n";
      for (TransformationDescription::DataPoints::const_iterator it = transformation.getDataPoints().begin();
           it != transformation.getDataPoints().end(); ++it)
      {
        os << "\t\t\t<Pair from=\"" << it->first << "\" to=\"" << it->second;
        if (!it->note.empty())
        {
          os << "\" note=\"" << writeXMLEscape(it->note);
        }
        os << "\"/>\n";
      }
      os << "\t\t</Pairs>\n";
    }

    // close tag
    os << "\t</Transformation>\n";

    //write footer
    os << "</TrafoXML>\n";

    //close stream
    os.close();
  }

  void TransformationXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {

    String element = sm_.convert(qname);

    if (element == "TrafoXML")
    {
      //check file version against schema version
      double file_version = attributeAsDouble_(attributes, "version");
      if (file_version > version_.toDouble())
      {
        warning(LOAD, String("The XML file (") + file_version + ") is newer than the parser (" + version_ + "). This might lead to undefined program behavior.");
      }
    }
    else if (element == "Transformation")
    {
      model_type_ = attributeAsString_(attributes, "name");
    }
    else if (element == "Param")
    {
      String type = attributeAsString_(attributes, "type");
      if (type == "int")
      {
        params_.setValue(attributeAsString_(attributes, "name"), attributeAsInt_(attributes, "value"));
      }
      else if (type == "float")
      {
        params_.setValue(attributeAsString_(attributes, "name"), attributeAsDouble_(attributes, "value"));
      }
      else if (type == "string")
      {
        params_.setValue(attributeAsString_(attributes, "name"), String(attributeAsString_(attributes, "value")));
      }
      else
      {
        error(LOAD, String("Unsupported parameter type: '") + type + "'");
      }

    }
    else if (element == "Pairs")
    {
      data_.reserve(attributeAsInt_(attributes, "count"));
    }
    else if (element == "Pair")
    {
      TransformationDescription::DataPoint point;
      point.first = attributeAsDouble_(attributes, "from");
      point.second = attributeAsDouble_(attributes, "to");
      optionalAttributeAsString_(point.note, attributes, "note");
      data_.push_back(point);
    }
    else
    {
      warning(LOAD, String("Unknown element: '") + element + "'");
    }
  }

} // namespace OpenMS
