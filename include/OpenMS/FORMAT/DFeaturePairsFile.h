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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_DFEATUREPAIRSFILE_H
#define OPENMS_FORMAT_DFEATUREPAIRSFILE_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/DFeaturePairsHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <iostream>  
#include <fstream>   
#include <sstream>
#include <vector>

namespace OpenMS
{
	
	/**
		@brief This class provides Input/Output functionality for the class DFeaturePairVector.
		
		The features pairs are computed by an instance of DBaseFeatureMatcher during the
		matching of MS maps. The features pairs are stored in a pseudo XML format. No schema
		has been developed yet therefore no validation can be performed.
  
  	@ingroup FileIO
  */
  class DFeaturePairsFile
  {
    public:
       /** @name Constructors and Destructor */
      //@{
      ///Default constructor
      DFeaturePairsFile() 
      {
      	
      }
       ///Destructor
      virtual ~DFeaturePairsFile() 
      {
      	
      }
      //@}

      /** @name Accessors */
      //@{
      /// loads the file with name @p filename into @p pairs.
      template<Size D> 
      void load(String filename, DFeaturePairVector<D>& pairs) throw (Exception::FileNotFound, Exception::ParseError)
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
				Internal::DFeaturePairsHandler<D> handler(pairs,filename);
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

      /// stores the pair vector @p pairs in file with name @p filename. 
      template<Size D> 
      void store(String filename, const DFeaturePairVector<D>& pairs) const throw (Exception::UnableToCreateFile)
      {
			  if (pairs.empty()) return;

				std::ofstream os(filename.c_str(),std::fstream::out);
				if (!os.is_open())
				{
					os.close();
					throw Exception::UnableToCreateFile(__FILE__,__LINE__,"DFeaturePairsFile::store()",filename);
				}

				Internal::DFeaturePairsHandler<D> handler(pairs,filename);
				handler.writeTo(os);
				os.close();
      }
               
      //@}

	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_DFEATUREPAIRSFILE_H
