// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Eva Lange$
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_CONSENSUSXMLFILE_H
#define OPENMS_FORMAT_CONSENSUSXMLFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/ConsensusXMLHandler.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>

namespace OpenMS
{
  /**
    @brief 

    @ingroup FileIO
  */
  class ConsensusXMLFile
  {
    public:
      ///Default constructor
      ConsensusXMLFile();
      ///Destructor
      ~ConsensusXMLFile();

      /**
          @brief Loads a consenus map from a ConsensusXML file.

          @p 
         */
      template <typename ElementT>
      void load(const String& filename, ConsensusMap<ElementT>& map) throw (Exception::FileNotFound, Exception::ParseError)
      {
        //try to open file
        if (!File::exists(filename))
        {
          throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
        }

        // initialize parser
        try
        {
          xercesc::XMLPlatformUtils::Initialize();
        }
        catch (const xercesc::XMLException& toCatch)
        {
          throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
        }

        xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
        parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
        parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);

        map = ConsensusMap<ElementT>();  // clear map
        Internal::ConsensusXMLHandler< StarAlignment<ElementT> > handler(map,filename);

        parser->setContentHandler(&handler);
        parser->setErrorHandler(&handler);

        xercesc::LocalFileInputSource source( xercesc::XMLString::transcode(filename.c_str()) );
        try
        {
          parser->parse(source);
          delete(parser);
        }
        catch (const xercesc::XMLException& toCatch)
        {
          throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
        }
        catch (const xercesc::SAXException& toCatch)
        {
          throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
        }
      }

      /**
        @brief Stores a staralignment object into consensusXML format.

        @p 
      */
      template <typename AlignmentT>
      void store(const String& filename, const AlignmentT& alignment)
      const throw (Exception::UnableToCreateFile)
      {
        std::ofstream os(filename.c_str());
        if (!os)
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__,__PRETTY_FUNCTION__,filename);
        }

        //read data and close stream
        Internal::ConsensusXMLHandler<AlignmentT> handler(alignment,filename);
        handler.writeTo(os);
        os.close();
      }
  };

} // namespace OpenMS

#endif // OPENMS_FOMAT_MZXMLFILE_H
