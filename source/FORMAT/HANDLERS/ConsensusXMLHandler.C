// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FORMAT/HANDLERS/ConsensusXMLHandler.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{
  namespace Internal
  {
    /// Load the peaks
    template <>
    template <>
    void ConsensusXMLHandler< StarAlignment< ConsensusFeature< FeatureMap > > >::loadFile_< ConsensusFeature< FeatureMap > >(const String& file_name, UnsignedInt id, const ConsensusFeature< FeatureMap >& /* c */ ) throw (Exception::FileNotFound, Exception::ParseError)
    {
      DFeatureMapFile feature_file;
      feature_file.load(file_name,(consensus_map_->getMapVector())[id]);
    }

    // load MzData
    template <>
    template <>
    void ConsensusXMLHandler< StarAlignment< ConsensusPeak< DPeakArray<2,Peak> > > >::loadFile_< ConsensusPeak< DPeakArray<2,Peak> > >( const String& file_name, UnsignedInt id, const ConsensusPeak< DPeakArray<2,Peak> >& /* c */) throw (Exception::FileNotFound, Exception::ParseError)
    {
      MzDataFile mzdata_file;
      MSExperiment< Peak > ms_exp;
      mzdata_file.load(file_name,ms_exp);
      ms_exp.get2DData((consensus_map_->getMapVector())[id]);
    }

    // load consensusXML
    template <>
    template <>
    void ConsensusXMLHandler< StarAlignment< ConsensusFeature< ConsensusMap< ConsensusFeature< FeatureMap > > > > >::loadFile_<ConsensusFeature< ConsensusMap< ConsensusFeature< FeatureMap > > > >(const String& file_name, UnsignedInt id, const ConsensusFeature< ConsensusMap< ConsensusFeature< FeatureMap > > >& /* c */) throw (Exception::FileNotFound, Exception::ParseError)
    {
      ConsensusXMLHandler< StarAlignment< ConsensusFeature<> > > handler(((consensus_map_->getMapVector())[id]),file_name);

      //try to open file
      if (!File::exists(file_name))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, file_name);
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

      parser->setContentHandler(&handler);
      parser->setErrorHandler(&handler);

      // try to parse file
      xercesc::LocalFileInputSource source( xercesc::XMLString::transcode(file_name.c_str()) );
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
  } // namespace Internal
} // namespace OpenMS






