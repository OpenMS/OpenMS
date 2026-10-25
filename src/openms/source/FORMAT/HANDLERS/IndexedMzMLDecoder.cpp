// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/CONCEPT/Init.h>

#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/DATAACCESS/S3InputSource.h>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>

#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/util/XMLString.hpp>

namespace OpenMS
{

  namespace IndexedMzMLUtils
  {
    std::streampos stringToStreampos(const std::string& s)
    {
      // Try to cast the string to a type that can hold the integer value
      // stored in the std::streampos type (which can vary from 32 to 128 bits
      // depending on the system).
      //
      // long long has a minimal size of 64 bits and can address a range up to
      // 16 Exbibit (or 2 Exbibyte), we can hopefully expect our files to be
      // smaller than an Exabyte and should be safe.
      std::streampos res;
      try
      {

        res = boost::lexical_cast< unsigned long long >(s);
        // TESTING CAST: res = (int) res; // use this only when emulating 32 bit systems

      }
      catch (boost::bad_lexical_cast&)
      {
        std::cerr << "Trying to convert corrupted / unreadable value to std::streampos : " << s << std::endl;
        std::cerr << "This can also happen if the value exceeds 63 bits, please check your input." << std::endl;
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            std::string("Could not convert string '") + s + "' to a 64 bit integer.");
      }

      // Check if the value can fit into std::streampos
      if ( std::abs( boost::lexical_cast< long double >(s) - res) > 0.1)
      {
        std::cerr << "Your system may not support addressing a file of this size,"
          << " only addresses that fit into a " << sizeof(std::streamsize)*8 <<
          " bit integer are supported on your system." << std::endl;
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            std::string("Could not convert string '") + s + "' to an integer on your system.");
      }

      return res;
    }
  }

  int IndexedMzMLDecoder::parseOffsets(const std::string& filename, std::streampos indexoffset, OffsetVector& spectra_offsets, OffsetVector& chromatograms_offsets)
  {
    std::string buffer_str;
    long long length = 0;

    if (S3InputSource::isS3URI(filename))
    {
      //-------------------------------------------------------------
      // S3 object: read from indexoffset to its end
      //-------------------------------------------------------------
      length = S3InputSource::objectSize(filename);

      if (indexoffset < 0 || indexoffset > length)
      {
        std::cerr << "IndexedMzMLDecoder::parseOffsets Error: Offset was " <<
          indexoffset << " (not between 0 and " << length << ")." << std::endl;
        return -1;
      }

      buffer_str = S3InputSource::readRange(filename, indexoffset, length - indexoffset);
    }
    else
    {
    //-------------------------------------------------------------
    // Open file, jump to end and read last indexoffset bytes into buffer.
    //-------------------------------------------------------------
    std::ifstream f(filename.c_str());
    if (!f.is_open())
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }

    // get length of file:
    f.seekg(0, f.end);
    length = f.tellg();

    if (indexoffset < 0 || indexoffset > length)
    {
      std::cerr << "IndexedMzMLDecoder::parseOffsets Error: Offset was " <<
        indexoffset << " (not between 0 and " << length << ")." << std::endl;
      return -1;
    }

    //-------------------------------------------------------------
    // Read full end of file to parse offsets for spectra and chroms
    //-------------------------------------------------------------
    // read data as a block into a buffer

    // allocate enough memory in buffer (+1 for string termination)
    std::streampos readl = length - indexoffset;
    char* buffer = new(std::nothrow) char[readl + std::streampos(1)];

    // catch case where not enough memory is available
    if (buffer == nullptr)
    {
      std::cerr << "IndexedMzMLDecoder::parseOffsets Could not allocate enough memory to read in index of indexedMzML" << std::endl;
      std::cerr << "IndexedMzMLDecoder::parseOffsets calculated index offset " << indexoffset << " and file length " << length <<
        ", consequently tried to read into memory " << readl << " bytes." << std::endl;
      return -1;
    }

    // read into memory
    f.seekg(-readl, f.end);
    f.read(buffer, readl);
    buffer[readl] = '\0';
    buffer_str = buffer;
    delete[] buffer;
    }

    //-------------------------------------------------------------
    // Add a sane start element and then give it to a DOM parser
    //-------------------------------------------------------------
    // http://stackoverflow.com/questions/4691039/making-xerces-parse-a-string-insted-of-a-file
    std::string tmp_fixed_xml = "<indexedmzML>" + buffer_str + "\n";
    int res = domParseIndexedEnd_(tmp_fixed_xml, spectra_offsets, chromatograms_offsets);

    return res;
  }

  std::streampos IndexedMzMLDecoder::findIndexListOffset(const std::string& filename, int buffersize)
  {
    // return value
    std::streampos indexoffset = -1;
    std::string buffer_str;

    if (S3InputSource::isS3URI(filename))
    {
      //-------------------------------------------------------------
      // S3 object: read its last buffersize bytes
      //-------------------------------------------------------------
      const Int64 size = S3InputSource::objectSize(filename);
      const Int64 start = std::max<Int64>(0, size - buffersize);
      buffer_str = S3InputSource::readRange(filename, start, size - start);
    }
    else
    {
    //-------------------------------------------------------------
    // Open file, jump to end and read last n bytes into buffer.
    //-------------------------------------------------------------
    std::ifstream f(filename.c_str());

    if (!f.is_open())
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }

    // Read the last few bytes and hope our offset is there to be found
    std::unique_ptr<char[]> buffer(new char[buffersize + 1]);
    f.seekg(-buffersize, f.end);
    f.read(buffer.get(), buffersize);
    buffer.get()[buffersize] = '\0';
    buffer_str = buffer.get();
    }

#ifdef DEBUG_READER
    std::cout << " reading file " << filename  << " with size " << buffersize << '\n';
    std::cout << buffer_str << '\n';
#endif

    //-------------------------------------------------------------
    // Since we could be anywhere in the XML structure, use regex to find
    // indexListOffset and read its content.
    //-------------------------------------------------------------
    boost::regex listoffset_rx(R"(<[^>/]*indexListOffset\s*>\s*(\d*))");
    boost::smatch matches;
    boost::regex_search(buffer_str, matches, listoffset_rx);
    std::string thismatch(matches[1].first, matches[1].second);
    if (!thismatch.empty())
    {
      try
      {
        indexoffset = OpenMS::IndexedMzMLUtils::stringToStreampos(thismatch);
      }
      catch (Exception::ConversionError& /*e*/)
      {
        std::cerr << "Corrupted / unreadable value in <indexListOffset> : " << thismatch << std::endl;
        // free resources and re-throw
        throw;  // re-throw conversion error
      }
    }
    else
    {
      std::cerr << "IndexedMzMLDecoder::findIndexListOffset Error: Could not find element indexListOffset in the last "
                << buffersize << " bytes. Maybe this is not a indexedMzML."
                << buffer_str << std::endl;
    }

    return indexoffset;
  }

  int IndexedMzMLDecoder::domParseIndexedEnd_(const std::string& in, OffsetVector& spectra_offsets, OffsetVector& chromatograms_offsets)
  {

    /*

     We parse something like 
     
      <indexedmzML>
        <indexList count="1">
          <index name="chromatogram">
            <offset idRef="1">9752</offset>
          </index>
        </indexList>
        <indexListOffset>26795</indexListOffset>
      <fileChecksum>0</fileChecksum>
      </indexedmzML>

    */

    //-------------------------------------------------------------
    // Create parser from input string using MemBufInputSource
    //-------------------------------------------------------------
    // Xerces has to be initialised; CONCEPT/Init.cpp does so on loading the library, but a
    // static library (WebAssembly) links that object file only if something else calls into it
    static const bool xerces_initialized = (Internal::xercesInitialize(), true);
    (void)xerces_initialized;
    xercesc::MemBufInputSource myxml_buf(
        reinterpret_cast<const unsigned char*>(in.c_str()), in.length(), "myxml (in memory)");
    xercesc::XercesDOMParser parser;
    parser.setDoNamespaces(false);
    parser.setDoSchema(false);
    parser.setLoadExternalDTD(false);
    parser.parse(myxml_buf);

    //-------------------------------------------------------------
    // Start parsing
    // see http://www.yolinux.com/TUTORIALS/XML-Xerces-C.html
    //-------------------------------------------------------------

    // no need to free this pointer - owned by the parent parser object
    xercesc::DOMDocument* doc =  parser.getDocument();
    // Get the top-level element ("indexedmzML")
    xercesc::DOMElement* elementRoot = doc->getDocumentElement();
    if (!elementRoot)
    {
      std::cerr << "IndexedMzMLDecoder::domParseIndexedEnd Error: " <<
        "No root element found:" << std::endl << std::endl << in << std::endl;
      return -1;
    }

    // Extract the indexList tag (there should only be one)
    XMLCh* x_tag = xercesc::XMLString::transcode("indexList");
    xercesc::DOMNodeList* li = elementRoot->getElementsByTagName(x_tag);
    xercesc::XMLString::release(&x_tag);
    if (li->getLength() != 1)
    {
      std::cerr << "IndexedMzMLDecoder::domParseIndexedEnd Error: "
        << "no indexList element found:" << std::endl << std::endl << in << std::endl;
      return -1;
    }
    xercesc::DOMNode* indexListNode = li->item(0);

    XMLCh* x_idref_tag = xercesc::XMLString::transcode("idRef");
    XMLCh* x_name_tag = xercesc::XMLString::transcode("name");

    xercesc::DOMNodeList* index_elems = indexListNode->getChildNodes();
    const  XMLSize_t nodeCount_ = index_elems->getLength();

    // Iterate through indexList elements (only two elements should be present
    // which should be either spectrum or chromatogram offsets)
    for (XMLSize_t j = 0; j < nodeCount_; ++j)
    {
      xercesc::DOMNode* currentNode = index_elems->item(j);
      if (currentNode->getNodeType() && // true is not NULL
          currentNode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE) // is element
      {
        std::vector<std::pair<std::string, std::streampos> > result;

        // Iterate through children
        // NOTE: Using xercesc::DOMNodeList and "item" is a very bad idea since
        //       each "item" call has complexity of O(n), see the
        //       implementation in DOMNodeListImpl.cpp :
        //       https://svn.apache.org/repos/asf/xerces/c/trunk/src/xercesc/dom/impl/DOMNodeListImpl.cpp
        //
        for (xercesc::DOMNode* currentONode = currentNode->getFirstChild(); currentONode != nullptr;
             currentONode = currentONode->getNextSibling())
        {
          if (currentONode->getNodeType() && // true is not NULL
              currentONode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE) // is element
          {
            xercesc::DOMElement* currentElement = dynamic_cast<xercesc::DOMElement*>(currentONode);

            char* x_name = xercesc::XMLString::transcode(currentElement->getAttribute(x_idref_tag));
            char* x_offset = xercesc::XMLString::transcode(currentONode->getTextContent());

            std::streampos thisOffset = OpenMS::IndexedMzMLUtils::stringToStreampos(std::string(x_offset) );
            result.emplace_back(std::string(x_name), thisOffset);

            xercesc::XMLString::release(&x_name);
            xercesc::XMLString::release(&x_offset);
          }
        }

        // should be either spectrum or chromatogram ...
        xercesc::DOMElement* currentElement = dynamic_cast<xercesc::DOMElement*>(currentNode);
        char* x_indexName = xercesc::XMLString::transcode(currentElement->getAttribute(x_name_tag));
        std::string name(x_indexName);
        xercesc::XMLString::release(&x_indexName);

        if (name == "spectrum")
        {
          spectra_offsets = result;
        }
        else if (name == "chromatogram")
        {
          chromatograms_offsets = result;
        }
        else
        {
          std::cerr << "IndexedMzMLDecoder::domParseIndexedEnd Error: expected only " <<
            "'spectrum' or 'chromatogram' below indexList but found instead '" <<
            name << "'." << std::endl;
          xercesc::XMLString::release(&x_idref_tag);
          xercesc::XMLString::release(&x_name_tag);
          return -1;
        }
      }
    }
    xercesc::XMLString::release(&x_idref_tag);
    xercesc::XMLString::release(&x_name_tag);

    return 0;
  }

}
