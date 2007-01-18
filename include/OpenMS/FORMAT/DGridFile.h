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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_DGRIDFILE_H
#define OPENMS_FORMAT_DGRIDFILE_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/DGridHandler.h>
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
  	@brief Provides Input/Output functionality for instances of class DGrid.
  	
  	@ingroup FileIO
  */
  class DGridFile
  {
    public:
       /** @name Constructors and Destructor */
      //@{
      ///Default constructor
      DGridFile() 
      { 
      	
      }
       ///Destructor
      virtual ~DGridFile() 
      { 
      	
      }
      //@}

      /** @name Accessors */
      //@{
      /// loads the file with name @p filename into @p grid. 
      template<Size D> 
      void load(String filename, DGrid<D>& grid) throw (Exception::FileNotFound,Exception::ParseError)
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
				Internal::DGridHandler<D> handler(grid,filename);
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

      /// stores the grid @p grid in file with name @p filename. 
      template<Size D> 
      void store(String filename, const DGrid<D>& grid) const throw (Exception::UnableToCreateFile)
      {
			  if (grid.empty()) return;

				std::ofstream os(filename.c_str(),std::fstream::out);
				if (!os.is_open())
				{
					os.close();
					throw Exception::UnableToCreateFile(__FILE__,__LINE__,"DGridFile::store()",filename);
				}

				Internal::DGridHandler<D> handler(grid,filename);
				handler.writeTo(os);
				os.close();
      }
               
      //@}
	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_DGRIDFILE_H
