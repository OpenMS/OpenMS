// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>
#include <OpenMS/CONCEPT/Init.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <libxml/parser.h>
#include <libxml/SAX2.h>
#include <memory>
#include <limits>

namespace OpenMS
{

  /// Small internal function to check the default data vectors
  void checkData_(std::vector<Internal::MzMLHandlerHelper::BinaryData>& data,
      SignedSize x_index, SignedSize int_index,
      bool x_precision_64, bool int_precision_64)
  {
    // Error if intensity or m/z (RT) is encoded as int32|64 - they should be float32|64!
    if ((!data[x_index].ints_32.empty()) || (!data[x_index].ints_64.empty()))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "", "Encoding m/z or RT array as integer is not allowed!");
    }
    if ((!data[int_index].ints_32.empty()) || (!data[int_index].ints_64.empty()))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "", "Encoding intensity array as integer is not allowed!");
    }

    Size mz_size = x_precision_64 ? data[x_index].floats_64.size() : data[x_index].floats_32.size();
    Size int_size = int_precision_64 ? data[int_index].floats_64.size() : data[int_index].floats_32.size();

    // Check if int-size and mz-size are equal
    if (mz_size != int_size)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "", "Error, intensity and m/z array length are unequal");
    }
  }

  inline void fillDataArray(const std::vector<Internal::MzMLHandlerHelper::BinaryData>& data,
                            const OpenMS::Interfaces::BinaryDataArrayPtr& array, bool precision_64, SignedSize index)
  {
    // This seems to be the fastest method to move the data (faster than copy or assign)
    if (precision_64)
    {
      array->data.insert(array->data.begin(), data[index].floats_64.begin(), data[index].floats_64.end());
    }
    else
    {
      array->data.insert(array->data.begin(), data[index].floats_32.begin(), data[index].floats_32.end());
    }
  }

  template<class T>
  inline void fillDataArrayFloat(const Internal::MzMLHandlerHelper::BinaryData& data, T& spectrum)
  {
    //create new array
    spectrum.getFloatDataArrays().resize(spectrum.getFloatDataArrays().size() + 1);
    //reserve space in the array
    spectrum.getFloatDataArrays().back().reserve(data.size);
    //copy meta info into MetaInfoDescription
    spectrum.getFloatDataArrays().back().MetaInfoDescription::operator=(data.meta);

    if (data.precision == Internal::MzMLHandlerHelper::BinaryData::PRE_64)
    {
      for (Size n = 0; n < data.floats_64.size(); n++)
      {
        double value = data.floats_64[n];
        spectrum.getFloatDataArrays().back().push_back(value);
      }
    }
    else
    {
      for (Size n = 0; n < data.floats_32.size(); n++)
      {
        double value = data.floats_32[n];
        spectrum.getFloatDataArrays().back().push_back(value);
      }
    }
  }

  template<class T>
  inline void fillDataArrayInt(const Internal::MzMLHandlerHelper::BinaryData& data, T& spectrum)
  {
    //create new array
    spectrum.getIntegerDataArrays().resize(spectrum.getIntegerDataArrays().size() + 1);
    //reserve space in the array
    spectrum.getIntegerDataArrays().back().reserve(data.size);
    //copy meta info into MetaInfoDescription
    spectrum.getIntegerDataArrays().back().MetaInfoDescription::operator=(data.meta);

    if (data.precision == Internal::MzMLHandlerHelper::BinaryData::PRE_64)
    {
      for (Size n = 0; n < data.ints_64.size(); n++)
      {
        double value = data.ints_64[n];
        spectrum.getIntegerDataArrays().back().push_back(value);
      }
    }
    else
    {
      for (Size n = 0; n < data.ints_32.size(); n++)
      {
        double value = data.ints_32[n];
        spectrum.getIntegerDataArrays().back().push_back(value);
      }
    }
  }

  template<class T>
  inline void fillDataArrayString(const Internal::MzMLHandlerHelper::BinaryData& data, T& spectrum)
  {
    //create new array
    spectrum.getStringDataArrays().resize(spectrum.getStringDataArrays().size() + 1);
    //reserve space in the array
    spectrum.getStringDataArrays().back().reserve(data.decoded_char.size());
    //copy meta info into MetaInfoDescription
    spectrum.getStringDataArrays().back().MetaInfoDescription::operator=(data.meta);

    if (data.precision == Internal::MzMLHandlerHelper::BinaryData::PRE_64)
    {
      for (Size n = 0; n < data.decoded_char.size(); n++)
      {
        std::string value = data.decoded_char[n];
        spectrum.getStringDataArrays().back().push_back(value);
      }
    }
  }

  inline void fillDataArrays(const std::vector<Internal::MzMLHandlerHelper::BinaryData>& input_data, MSSpectrum& spectrum)
  {
    for (Size i = 0; i < input_data.size(); i++) //loop over all binary data arrays
    {
      if (input_data[i].meta.getName() != "m/z array" && input_data[i].meta.getName() != "intensity array") // is meta data array?
      {
        if (input_data[i].data_type == Internal::MzMLHandlerHelper::BinaryData::DT_FLOAT)
        {
          fillDataArrayFloat(input_data[i], spectrum);
        }
        else if (input_data[i].data_type == Internal::MzMLHandlerHelper::BinaryData::DT_INT)
        {
          fillDataArrayInt(input_data[i], spectrum);
        }
        else if (input_data[i].data_type == Internal::MzMLHandlerHelper::BinaryData::DT_STRING)
        {
          fillDataArrayString(input_data[i], spectrum);
        }
      }
    }
  }

  inline void fillDataArrays(const std::vector<Internal::MzMLHandlerHelper::BinaryData>& input_data, MSChromatogram& spectrum)
  {
    for (Size i = 0; i < input_data.size(); i++) //loop over all binary data arrays
    {
      if (input_data[i].meta.getName() != "time array" && input_data[i].meta.getName() != "intensity array") // is meta data array?
      {
        if (input_data[i].data_type == Internal::MzMLHandlerHelper::BinaryData::DT_FLOAT)
        {
          fillDataArrayFloat(input_data[i], spectrum);
        }
        else if (input_data[i].data_type == Internal::MzMLHandlerHelper::BinaryData::DT_INT)
        {
          fillDataArrayInt(input_data[i], spectrum);
        }
        else if (input_data[i].data_type == Internal::MzMLHandlerHelper::BinaryData::DT_STRING)
        {
          fillDataArrayString(input_data[i], spectrum);
        }
      }
    }
  }

  template<class T>
  inline void fillDataArray(const std::vector<Internal::MzMLHandlerHelper::BinaryData>& data,
                            T& peak_container,
                            const bool x_precision_64,
                            const bool int_precision_64,
                            const SignedSize x_index,
                            const SignedSize int_index,
                            const Size default_array_length_)
  {
    typename T::PeakType tmp;
    if (x_precision_64 && !int_precision_64)
    {
      std::vector< double >::const_iterator mz_it = data[x_index].floats_64.begin();
      std::vector< float >::const_iterator int_it = data[int_index].floats_32.begin();
      for (Size n = 0; n < default_array_length_; n++)
      {
        //add peak
        tmp.setIntensity(*int_it);
        tmp.setPos(*mz_it);
        ++mz_it;
        ++int_it;
        peak_container.push_back(tmp);
      }
    }
    else if (x_precision_64 && int_precision_64)
    {
      std::vector< double >::const_iterator mz_it = data[x_index].floats_64.begin();
      std::vector< double >::const_iterator int_it = data[int_index].floats_64.begin();
      for (Size n = 0; n < default_array_length_; n++)
      {
        //add peak
        tmp.setIntensity(*int_it);
        tmp.setPos(*mz_it);
        ++mz_it;
        ++int_it;
        peak_container.push_back(tmp);
      }
    }
    else if (!x_precision_64 && int_precision_64)
    {
      std::vector< float >::const_iterator mz_it = data[x_index].floats_32.begin();
      std::vector< double >::const_iterator int_it = data[int_index].floats_64.begin();
      for (Size n = 0; n < default_array_length_; n++)
      {
        //add peak
        tmp.setIntensity(*int_it);
        tmp.setPos(*mz_it);
        ++mz_it;
        ++int_it;
        peak_container.push_back(tmp);
      }
    }
    else if (!x_precision_64 && !int_precision_64)
    {
      std::vector< float >::const_iterator mz_it = data[x_index].floats_32.begin();
      std::vector< float >::const_iterator int_it = data[int_index].floats_32.begin();
      for (Size n = 0; n < default_array_length_; n++)
      {
        //add peak
        tmp.setIntensity(*int_it);
        tmp.setPos(*mz_it);
        ++mz_it;
        ++int_it;
        peak_container.push_back(tmp);
      }
    }
  }

  void MzMLSpectrumDecoder::decodeBinaryDataMSSpectrum_(std::vector<BinaryData>& data, OpenMS::MSSpectrum& spectrum) const
  {
    Internal::MzMLHandlerHelper::decodeBase64Arrays(data, skip_xml_checks_);

    //look up the precision and the index of the intensity and m/z array
    bool x_precision_64 = true;
    bool int_precision_64 = true;
    SignedSize x_index = -1;
    SignedSize int_index = -1;
    Internal::MzMLHandlerHelper::computeDataProperties_(data, x_precision_64, x_index, "m/z array");
    Internal::MzMLHandlerHelper::computeDataProperties_(data, int_precision_64, int_index, "intensity array");

    //Abort if no m/z or intensity array is present
    if (int_index == -1 || x_index == -1)
    {
      std::cerr << "Error, intensity or m/z array is missing, skipping this spectrum" << std::endl;
      return;
    }

    checkData_(data, x_index, int_index, x_precision_64, int_precision_64);

    // At this point, we could check whether the defaultArrayLength from the
    // spectrum/chromatogram tag is equivalent to size of the decoded data
    Size default_array_length_ = x_precision_64 ? data[x_index].floats_64.size() : data[x_index].floats_32.size();

    spectrum.reserve(default_array_length_);
    fillDataArray(data, spectrum,
                  x_precision_64, int_precision_64,
                  x_index, int_index, default_array_length_);

    if (data.size() > 2) 
    {
      fillDataArrays(data, spectrum);
    }
  }

  void MzMLSpectrumDecoder::decodeBinaryDataMSChrom_(std::vector<BinaryData>& data, OpenMS::MSChromatogram& chromatogram) const
  {
    Internal::MzMLHandlerHelper::decodeBase64Arrays(data, skip_xml_checks_);

    //look up the precision and the index of the intensity and m/z array
    bool x_precision_64 = true;
    bool int_precision_64 = true;
    SignedSize x_index = -1;
    SignedSize int_index = -1;
    Internal::MzMLHandlerHelper::computeDataProperties_(data, x_precision_64, x_index, "time array");
    Internal::MzMLHandlerHelper::computeDataProperties_(data, int_precision_64, int_index, "intensity array");

    //Abort if no m/z or intensity array is present
    if (int_index == -1 || x_index == -1)
    {
      std::cerr << "Error, intensity or RT array is missing, skipping this spectrum" << std::endl;
      return;
    }

    checkData_(data, x_index, int_index, x_precision_64, int_precision_64);

    // At this point, we could check whether the defaultArrayLength from the
    // spectrum/chromatogram tag is equivalent to size of the decoded data
    Size default_array_length_ = x_precision_64 ? data[x_index].floats_64.size() : data[x_index].floats_32.size();

    chromatogram.reserve(default_array_length_);
    fillDataArray(data, chromatogram,
                  x_precision_64, int_precision_64,
                  x_index, int_index, default_array_length_);

    if (data.size() > 2) 
    {
      fillDataArrays(data, chromatogram);
    }
  }

  OpenMS::Interfaces::SpectrumPtr MzMLSpectrumDecoder::decodeBinaryDataSpectrum_(std::vector<BinaryData>& data) const
  {
    Internal::MzMLHandlerHelper::decodeBase64Arrays(data, skip_xml_checks_);
    OpenMS::Interfaces::SpectrumPtr sptr(new OpenMS::Interfaces::Spectrum);

    //look up the precision and the index of the intensity and m/z array
    bool x_precision_64 = true;
    bool int_precision_64 = true;
    SignedSize x_index = -1;
    SignedSize int_index = -1;
    Internal::MzMLHandlerHelper::computeDataProperties_(data, x_precision_64, x_index, "m/z array");
    Internal::MzMLHandlerHelper::computeDataProperties_(data, int_precision_64, int_index, "intensity array");

    //Abort if no m/z or intensity array is present
    if (int_index == -1 || x_index == -1)
    {
      std::cerr << "Error, intensity or m/z array is missing, skipping this spectrum" << std::endl;
      return sptr;
    }

    checkData_(data, x_index, int_index, x_precision_64, int_precision_64);

    // At this point, we could check whether the defaultArrayLength from the
    // spectrum/chromatogram tag is equivalent to size of the decoded data
    Size default_array_length_ = x_precision_64 ? data[x_index].floats_64.size() : data[x_index].floats_32.size();

    // TODO: also handle non-default arrays
    if (data.size() > 2)
    {
      std::cout << "MzMLSpectrumDecoder currently cannot handle meta data arrays, they are ignored." << std::endl;
    }

    // TODO: handle meta data from the binaryDataArray tag -> currently ignored
    // since we have no place to store them
    // TODO: would need to adopt OpenMS::Interfaces::SpectrumPtr to store additional arrays
    OpenMS::Interfaces::BinaryDataArrayPtr intensity_array(new OpenMS::Interfaces::BinaryDataArray);
    OpenMS::Interfaces::BinaryDataArrayPtr x_array(new OpenMS::Interfaces::BinaryDataArray);
    x_array->data.reserve(default_array_length_);
    intensity_array->data.reserve(default_array_length_);

    fillDataArray(data, x_array, x_precision_64, x_index);
    fillDataArray(data, intensity_array, int_precision_64, int_index);

    // TODO the other arrays

    sptr->setMZArray(x_array);
    sptr->setIntensityArray(intensity_array);
    return sptr;
  }

  OpenMS::Interfaces::ChromatogramPtr MzMLSpectrumDecoder::decodeBinaryDataChrom_(std::vector<BinaryData>& data) const
  {
    Internal::MzMLHandlerHelper::decodeBase64Arrays(data, skip_xml_checks_);
    OpenMS::Interfaces::ChromatogramPtr sptr(new OpenMS::Interfaces::Chromatogram);

    //look up the precision and the index of the intensity and m/z array
    bool x_precision_64 = true;
    bool int_precision_64 = true;
    SignedSize x_index = -1;
    SignedSize int_index = -1;
    Internal::MzMLHandlerHelper::computeDataProperties_(data, x_precision_64, x_index, "time array");
    Internal::MzMLHandlerHelper::computeDataProperties_(data, int_precision_64, int_index, "intensity array");

    //Abort if no m/z or intensity array is present
    if (int_index == -1 || x_index == -1)
    {
      std::cerr << "Error, intensity or RT array is missing, skipping this spectrum" << std::endl;
      return sptr;
    }

    checkData_(data, x_index, int_index, x_precision_64, int_precision_64);

    // At this point, we could check whether the defaultArrayLength from the
    // spectrum/chromatogram tag is equivalent to size of the decoded data
    Size default_array_length_ = x_precision_64 ? data[x_index].floats_64.size() : data[x_index].floats_32.size();

    // TODO: also handle non-default arrays
    if (data.size() > 2)
    {
      std::cout << "MzMLSpectrumDecoder currently cannot handle meta data arrays, they are ignored." << std::endl;
    }

    // TODO: handle meta data from the binaryDataArray tag -> currently ignored
    // since we have no place to store them

    OpenMS::Interfaces::BinaryDataArrayPtr intensity_array(new OpenMS::Interfaces::BinaryDataArray);
    OpenMS::Interfaces::BinaryDataArrayPtr x_array(new OpenMS::Interfaces::BinaryDataArray);
    x_array->data.reserve(default_array_length_);
    intensity_array->data.reserve(default_array_length_);

    fillDataArray(data, x_array, x_precision_64, x_index);
    fillDataArray(data, intensity_array, int_precision_64, int_index);

    // TODO the other arrays

    sptr->setTimeArray(x_array);
    sptr->setIntensityArray(intensity_array);
    return sptr;
  }

  namespace
  {
    std::string xmlProperty(xmlNode* node, const char* name)
    {
      xmlChar* value = xmlGetProp(node, BAD_CAST name);
      const std::string result = value ? reinterpret_cast<const char*>(value) : "";
      xmlFree(value);
      return result;
    }
  }

  void MzMLSpectrumDecoder::handleBinaryDataArray_(void* opaque, std::vector<BinaryData>& data)
  {
    auto* array = static_cast<xmlNode*>(opaque);
    data.emplace_back();
    bool has_binary = false;
    for (xmlNode* node = array->children; node; node = node->next)
    {
      if (node->type != XML_ELEMENT_NODE) { continue; }
      if (xmlStrEqual(node->name, BAD_CAST "binary"))
      {
        has_binary = true;
        if (!node->children) { continue; }
        if (node->children->next || node->children->type != XML_TEXT_NODE)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "", "Invalid XML: 'binary' element can only have a single, text node child element.");
        }
        data.back().base64.append(reinterpret_cast<const char*>(node->children->content));
      }
      else if (xmlStrEqual(node->name, BAD_CAST "cvParam"))
      {
        Internal::MzMLHandlerHelper::handleBinaryDataArrayCVParam(data, xmlProperty(node, "accession"),
          xmlProperty(node, "value"), xmlProperty(node, "name"), xmlProperty(node, "unitAccession"));
      }
      else if (xmlStrEqual(node->name, BAD_CAST "userParam")) { std::cout << " unhandled userParam\n"; }
      else if (xmlStrEqual(node->name, BAD_CAST "referenceableParamGroupRef")) { std::cout << " unhandled referenceableParamGroupRef\n"; }
    }
    if (!has_binary)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Invalid XML: 'binary' element needs to be present at least once inside 'binaryDataArray' element.");
    }
  }

  std::string MzMLSpectrumDecoder::domParseString_(const std::string& in, std::vector<BinaryData>& data)
  {
    if (in.size() > static_cast<Size>(std::numeric_limits<int>::max()))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "XML fragment exceeds libxml2 memory API size limit");
    }
    std::unique_ptr<xmlParserCtxt, decltype(&xmlFreeParserCtxt)> parser(xmlNewParserCtxt(), xmlFreeParserCtxt);
    if (!parser) { throw std::bad_alloc(); }
    // Large spectra can have a base64 text node exceeding the default 10 MB cap.
    // HUGE also raises depth/name caps; external entities stay disabled and the
    // independent entity-amplification limit is explicitly retained at 5.
    xmlCtxtSetMaxAmplification(parser.get(), 5);
    // Indexed offset spans include the completed root followed by parent closing
    // tags (and sometimes the next list opening). Stop at the real XML root end,
    // using the parser rather than a lexical search through comments or CDATA.
    xmlDoc* completed_document = nullptr;
    parser->_private = &completed_document;
    parser->sax->endElementNs = [](void* opaque, const xmlChar* name, const xmlChar* prefix, const xmlChar* uri)
    {
      auto* ctxt = static_cast<xmlParserCtxt*>(opaque);
      xmlSAX2EndElementNs(ctxt, name, prefix, uri);
      if (!ctxt->node)
      {
        *static_cast<xmlDoc**>(ctxt->_private) = ctxt->myDoc;
        ctxt->myDoc = nullptr; // transfer ownership before xmlCtxtReadMemory handles USER_STOP
        xmlStopParser(ctxt);
      }
    };
    xmlDoc* returned_document = xmlCtxtReadMemory(parser.get(), in.data(), static_cast<int>(in.size()), "spectrum.xml", nullptr,
      XML_PARSE_NONET | XML_PARSE_NO_XXE | XML_PARSE_NOENT | XML_PARSE_HUGE | XML_PARSE_NOERROR | XML_PARSE_NOWARNING);
    std::unique_ptr<xmlDoc, decltype(&xmlFreeDoc)> document(completed_document ? completed_document : returned_document, xmlFreeDoc);
    xmlNode* root = document ? xmlDocGetRootElement(document.get()) : nullptr;
    if (!root) { throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, in, "No root element or malformed XML"); }
    OPENMS_PRECONDITION(xmlStrEqual(root->name, BAD_CAST "spectrum") || xmlStrEqual(root->name, BAD_CAST "chromatogram"),
      "The input needs to contain a <spectrum> or <chromatogram> tag as root element.")
    if (!xmlHasProp(root, BAD_CAST "defaultArrayLength"))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, in, "Root element does not contain defaultArrayLength XML tag.");
    }
    const int length = StringUtils::toInt32(xmlProperty(root, "defaultArrayLength"));
    // Visit descendant binary arrays without allocating a second node list.
    for (xmlNode* node = root->children; node;)
    {
      if (node->type == XML_ELEMENT_NODE && xmlStrEqual(node->name, BAD_CAST "binaryDataArray"))
      {
        handleBinaryDataArray_(node, data);
        data.back().size = length;
      }
      if (node->children) { node = node->children; continue; }
      while (node != root && !node->next) { node = node->parent; }
      if (node == root) { break; }
      node = node->next;
    }
    return xmlProperty(root, "id");
  }

  void MzMLSpectrumDecoder::domParseSpectrum(const std::string& in, OpenMS::Interfaces::SpectrumPtr& sptr)
  {
    std::vector<BinaryData> data;
    domParseString_(in, data);
    sptr = decodeBinaryDataSpectrum_(data);
  }

  void MzMLSpectrumDecoder::domParseSpectrum(const std::string& in, MSSpectrum& s)
  {
    std::vector<BinaryData> data;
    std::string id = domParseString_(in, data);
    decodeBinaryDataMSSpectrum_(data, s);
    s.setNativeID(id);
  }

  void MzMLSpectrumDecoder::domParseChromatogram(const std::string& in, MSChromatogram& c)
  {
    std::vector<BinaryData> data;
    std::string id = domParseString_(in, data);
    decodeBinaryDataMSChrom_(data, c);
    c.setNativeID(id);
  }

  void MzMLSpectrumDecoder::domParseChromatogram(const std::string& in, OpenMS::Interfaces::ChromatogramPtr& sptr)
  {
    std::vector<BinaryData> data;
    domParseString_(in, data);
    sptr = decodeBinaryDataChrom_(data);
  }

  void MzMLSpectrumDecoder::setSkipXMLChecks(bool skip)
  {
    skip_xml_checks_ = skip;
  }

}

