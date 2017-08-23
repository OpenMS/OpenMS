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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
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

  void TransformationXMLFile::store(String filename, const TransformationDescription& transformation)
  {
    if (transformation.getModelType() == "")
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
    os << "<TrafoXML version=\"" << getVersion() << "\" xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/"
       << schema_location_.suffix('/') << "\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

    // open tag
    os << "\t<Transformation name=\"" << transformation.getModelType()
       << "\">\n";

    // write parameters
    Param params = transformation.getModelParameters();
    for (Param::ParamIterator it = params.begin(); it != params.end(); ++it)
    {
      if (it->value.valueType() != DataValue::EMPTY_VALUE)
      {
        switch (it->value.valueType())
        {
        case DataValue::INT_VALUE:
          os << "\t\t<Param  type=\"int\" name=\"" << it->name << "\" value=\"" << it->value.toString() << "\"/>\n";
          break;

        case DataValue::DOUBLE_VALUE:
          os << "\t\t<Param  type=\"float\" name=\"" << it->name << "\" value=\"" << it->value.toString() << "\"/>\n";
          break;

        case DataValue::STRING_VALUE:
        case DataValue::STRING_LIST:
        case DataValue::INT_LIST:
        case DataValue::DOUBLE_LIST:
          os << "\t\t<Param  type=\"string\" name=\"" << it->name << "\" value=\"" << it->value.toString() << "\"/>\n";
          break;

        default:         // no other value types are supported!
          fatalError(STORE, String("Unsupported parameter type of parameter '") + it->name + "' with value '" + it->value.toString() + "'");
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
