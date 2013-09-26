// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>

#include <boost/regex.hpp>
#include <fstream>
#include <string>

#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/util/XMLString.hpp>

namespace OpenMS
{
  int IndexedMzMLDecoder::parseOffsets(String in, int indexoffset, OffsetVector& spectra_offsets, OffsetVector& chromatograms_offsets)
  {
    //-------------------------------------------------------------
    // Open file, jump to end and read last indexoffset bytes into buffer.
    //-------------------------------------------------------------
    std::ifstream f(in.c_str());
    // get length of file:
    f.seekg(0, f.end);
    int length = f.tellg();

    if (indexoffset < 0 || indexoffset > length)
    {
      return -1;
    }

    //-------------------------------------------------------------
    // Read full end of file to parse offsets for spectra and chroms
    //-------------------------------------------------------------
    // read data as a block:
    // allocate memory:
    int readl = length - indexoffset;
    char* buffer = new char[readl + 1];
    f.seekg(-readl, f.end);
    f.read(buffer, readl);
    buffer[readl] = '\0';

    //-------------------------------------------------------------
    // Add a sane start element and then give it to a DOM parser
    //-------------------------------------------------------------
    // http://stackoverflow.com/questions/4691039/making-xerces-parse-a-string-insted-of-a-file
    std::string tmp_fixed_xml = "<indexedmzML>" +  String(buffer) + "\n";
    int res = domParseIndexedEnd(tmp_fixed_xml, spectra_offsets, chromatograms_offsets);

    delete[] buffer;

    return res;
  }

  int IndexedMzMLDecoder::findIndexListOffset(String in, int buffersize)
  {
    // return value
    int indexoffset = -1;

    //-------------------------------------------------------------
    // Open file, jump to end and read last n bytes into buffer.
    //-------------------------------------------------------------
    std::ifstream f(in.c_str());

    if (!f.is_open())
    {
      return indexoffset;
    }

    // Read the last few bytes and hope our offset is there to be found
    char* buffer = new char[buffersize + 1];
    f.seekg(-buffersize, f.end);
    f.read(buffer, buffersize);
    buffer[buffersize] = '\0';

#ifdef DEBUG_READER
    std::cout << " reading file " << in  << " with size " << buffersize << std::endl;
    std::cout << buffer << std::endl;
#endif

    //-------------------------------------------------------------
    // Since we could be anywhere in the XML structure, use regex to find
    // indexListOffset and read its content.
    //-------------------------------------------------------------
    std::string text(buffer);
    boost::regex listoffset_rx("<[^>/]*indexListOffset\\s*>\\s*(\\d*)");
    boost::cmatch matches;
    boost::regex_search(buffer, matches, listoffset_rx);
    String thismatch(matches[1].first, matches[1].second);
    if (thismatch.size() > 0)
    {
      // TODO catch if it throws
      indexoffset = thismatch.toInt();
    }
    else
    {
      std::cout << "Something is wrong here, could not find indexListOffset -- maybe this is not a indexedmzML? " << std::endl;
    }

    f.close();
    delete[] buffer;

    return indexoffset;
  }

  int IndexedMzMLDecoder::domParseIndexedEnd(std::string in, OffsetVector& spectra_offsets, OffsetVector& chromatograms_offsets)
  {
    // see http://www.yolinux.com/TUTORIALS/XML-Xerces-C.html
    xercesc::MemBufInputSource myxml_buf(reinterpret_cast<const unsigned char*>(in.c_str()), in.length(), "myxml (in memory)");
    xercesc::XercesDOMParser* parser = new xercesc::XercesDOMParser();
    parser->setDoNamespaces(false);
    parser->setDoSchema(false);
    parser->setLoadExternalDTD(false);
    parser->parse(myxml_buf);

    // no need to free this pointer - owned by the parent parser object
    xercesc::DOMDocument* doc =  parser->getDocument();
    // Get the top-level element ("indexedmzML")
    xercesc::DOMElement* elementRoot = doc->getDocumentElement();
    if (!elementRoot)
    {
      delete parser;
      return -1;
    }

    // Extract the indexList tag (there should only be one)
    XMLCh* tag = xercesc::XMLString::transcode("indexList");
    xercesc::DOMNodeList* li = elementRoot->getElementsByTagName(tag);
    xercesc::XMLString::release(&tag);
    if (li->getLength() != 1)
    {
      delete parser;
      return -1;
    }
    xercesc::DOMNode* indexListNode = li->item(0);

    // Iterate through indexList elements
    xercesc::DOMNodeList* index_elems = indexListNode->getChildNodes();
    const  XMLSize_t nodeCount_ = index_elems->getLength();
    for (XMLSize_t j = 0; j < nodeCount_; ++j)
    {
      xercesc::DOMNode* currentNode = index_elems->item(j);
      if (currentNode->getNodeType() && // true is not NULL
          currentNode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE) // is element
      {
        std::vector<std::pair<std::string, long> > result;
        xercesc::DOMNodeList* offset_elems = currentNode->getChildNodes();
        for (XMLSize_t k = 0; k < offset_elems->getLength(); ++k)
        {
          xercesc::DOMNode* currentONode = offset_elems->item(k);
          if (currentONode->getNodeType() && // true is not NULL
              currentONode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE) // is element
          {
            xercesc::DOMElement* currentElement = dynamic_cast<xercesc::DOMElement*>(currentONode);
            std::string name = xercesc::XMLString::transcode(currentElement->getAttribute(xercesc::XMLString::transcode("idRef")));
            long thisOffset = String(xercesc::XMLString::transcode(currentONode->getTextContent())).toInt();
            result.push_back(std::make_pair(name, thisOffset));
          }
        }

        // should be either spectrum or chromatogram ...
        xercesc::DOMElement* currentElement = dynamic_cast<xercesc::DOMElement*>(currentNode);
        std::string name = xercesc::XMLString::transcode(currentElement->getAttribute(xercesc::XMLString::transcode("name")));
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
          // TODO
          //throw 1;
          delete parser;
          return -1;
        }
      }
    }

    delete parser;
    return 0;
  }

}
