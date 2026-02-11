// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOMText.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include <memory>

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
        String value = data.decoded_char[n];
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

  void MzMLSpectrumDecoder::handleBinaryDataArray_(xercesc::DOMNode* indexListNode, std::vector<BinaryData>& data)
  {
    // access result through data.back()
    data.emplace_back();

    // using CONST_XMLCH since writing a u16 char array each time is ugly. disadvantage = no constexpr.
    // Uses only reinterpret_cast, so usually no runtime cost though.
    static const XMLCh* TAG_CV = CONST_XMLCH("cvParam");
    static const XMLCh* TAG_binary = CONST_XMLCH("binary");
    static const XMLCh* TAG_userParam = CONST_XMLCH("userParam");
    static const XMLCh* TAG_referenceableParamGroupRef = CONST_XMLCH("referenceableParamGroupRef");
    static const XMLCh* TAG_accession = CONST_XMLCH("accession");
    static const XMLCh* TAG_unit_accession = CONST_XMLCH("unitAccession");
    static const XMLCh* TAG_value = CONST_XMLCH("value");
    static const XMLCh* TAG_name = CONST_XMLCH("name");

    OpenMS::Internal::StringManager sm;

    // Iterate through binaryDataArray elements
    // only allowed subelements:
    //  - referenceableParamGroupRef (0+)
    //  - cvParam (0+)
    //  - userParam (0+)
    //  - binary (1)
    xercesc::DOMNodeList* index_elems = indexListNode->getChildNodes();
    const  XMLSize_t nodeCount_ = index_elems->getLength();
    bool has_binary_tag = false;
    for (XMLSize_t j = 0; j < nodeCount_; ++j)
    {
      xercesc::DOMNode* currentNode = index_elems->item(j);
      if (currentNode->getNodeType() &&   // true is not NULL
          currentNode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE)   // is element
      {
        xercesc::DOMElement* currentElement = dynamic_cast<xercesc::DOMElement*>(currentNode);
        if (xercesc::XMLString::equals(currentElement->getTagName(), TAG_binary))
        {
          // Found the <binary> tag
          has_binary_tag = true;

          // Skip any empty <binary></binary> tags
          if (!currentNode->hasChildNodes())
          {
            continue;
          }

          // Valid mzML does not have any other child nodes except text
          if (currentNode->getChildNodes()->getLength() != 1)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "", "Invalid XML: 'binary' element can only have a single, text node child element.");
          }

          // Now we know that the <binary> node has exactly one single child node, the text that we want!
          xercesc::DOMNode* textNode_ = currentNode->getFirstChild();

          if (textNode_->getNodeType() == xercesc::DOMNode::TEXT_NODE)
          {
            xercesc::DOMText* textNode (static_cast<xercesc::DOMText*> (textNode_));
            sm.appendASCII(textNode->getData(),
                textNode->getLength(), data.back().base64);
          }
          else
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "", "Invalid XML: 'binary' element can only have a single, text node child element.");
          }
        }
        else if (xercesc::XMLString::equals(currentElement->getTagName(), TAG_CV))
        {
          std::string accession = sm.convert(currentElement->getAttribute(TAG_accession));
          std::string value = sm.convert(currentElement->getAttribute(TAG_value));
          std::string name = sm.convert(currentElement->getAttribute(TAG_name));
          std::string unit_accession = sm.convert(currentElement->getAttribute(TAG_unit_accession));

          // set precision, data_type
          Internal::MzMLHandlerHelper::handleBinaryDataArrayCVParam(data, accession, value, name, unit_accession);
        }
        else if (xercesc::XMLString::equals(currentElement->getTagName(), TAG_userParam))
        {
          std::cout << " unhandled userParam" << std::endl;
        }
        else if (xercesc::XMLString::equals(currentElement->getTagName(), TAG_referenceableParamGroupRef))
        {
          std::cout << " unhandled referenceableParamGroupRef" << std::endl;
        }
        else
        {
          //std::cout << "unhandled" << (string)xercesc::XMLString::transcode(currentNode->getNodeName() << std::endl;
        }
      }
    }

    // Throw exception upon invalid mzML: the <binary> tag is required inside <binaryDataArray>
    if (!has_binary_tag)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "", "Invalid XML: 'binary' element needs to be present at least once inside 'binaryDataArray' element.");
    }
  }

  std::string MzMLSpectrumDecoder::domParseString_(const std::string& in, std::vector<BinaryData>& data)
  {
    // PRECONDITON is below (since we first need to do XML parsing before validating)
    // initializer list of XMLCh (= usually some type that fits utf16) from ASCII chars
    static constexpr XMLCh id_tag[] = {'i','d', 0};
    static constexpr XMLCh default_array_length_tag[] = { 'd','e','f','a','u','l','t','A','r','r','a','y','L','e','n','g','t','h', 0};
    static constexpr XMLCh binary_data_array_tag[] = { 'b','i','n','a','r','y','D','a','t','a','A','r','r','a','y', 0};

    //-------------------------------------------------------------
    // Create parser from input string using MemBufInputSource
    //-------------------------------------------------------------
    xercesc::MemBufInputSource myxml_buf(reinterpret_cast<const unsigned char*>(in.c_str()), in.length(), "myxml (in memory)");
    std::unique_ptr<xercesc::XercesDOMParser> parser = std::make_unique<xercesc::XercesDOMParser>(); // make sure it is deleted even in case of exceptions
    parser->setDoNamespaces(false);
    parser->setDoSchema(false);
    parser->setLoadExternalDTD(false);
    parser->parse(myxml_buf);

    //-------------------------------------------------------------
    // Start parsing
    // see http://www.yolinux.com/TUTORIALS/XML-Xerces-C.html
    //-------------------------------------------------------------

    // no need to free this pointer - owned by the parent parser object
    xercesc::DOMDocument* doc =  parser->getDocument();
    // Get the top-level element (needs to be <spectrum> or <chromatogram>)
    xercesc::DOMElement* elementRoot = doc->getDocumentElement();
    if (!elementRoot)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, in, "No root element");
    }

    OPENMS_PRECONDITION(xercesc::XMLString::equals(elementRoot->getTagName(), CONST_XMLCH("spectrum")) || xercesc::XMLString::equals(elementRoot->getTagName(), CONST_XMLCH("chromatogram")),
          (String("The input needs to contain a <spectrum> or <chromatogram> tag as root element. Got instead '") +
          String(Internal::unique_xerces_ptr(xercesc::XMLString::transcode(elementRoot->getTagName())).get()) + String("'.")).c_str())

    // defaultArrayLength is a required attribute for the spectrum and the
    // chromatogram tag (but still check for it first to be safe).
    if (elementRoot->getAttributeNode(default_array_length_tag) == nullptr)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          in, "Root element does not contain defaultArrayLength XML tag.");
    }
    int default_array_length = xercesc::XMLString::parseInt(elementRoot->getAttribute(default_array_length_tag));
    OpenMS::Internal::StringManager sm;
    std::string id = sm.convert(elementRoot->getAttribute(id_tag));

    // Extract the binaryDataArray elements (there may be multiple) and process them
    xercesc::DOMNodeList* li = elementRoot->getElementsByTagName(binary_data_array_tag);
    for (Size i = 0; i < li->getLength(); i++)
    {
      // Will append one single BinaryData object to data      
      handleBinaryDataArray_(li->item(i), data);
      
      // Set the size correctly (otherwise MzMLHandlerHelper complains).
      data.back().size = default_array_length;
    }

    return id;
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

