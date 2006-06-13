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
// $Id: MzXMLFile.h,v 1.21 2006/05/23 15:37:10 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MZXMLFILE_H
#define OPENMS_FORMAT_MZXMLFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/MzXMLHandler.h>

#include <fstream>

#include <qtextstream.h>
#include <qstring.h>

namespace OpenMS
{
	class String;
  /**
  	@brief File adapter for MzXML files
  
  	
  
  	@ingroup FileIO
  */
  class MzXMLFile
  {
		public:
      ///Default constructor
      MzXMLFile();
      ///Destructor
      ~MzXMLFile();


			/**
      	@brief Loads a map from a MzXML file.

      	@p map has to be a MSExperiment or have the same interface.
      */
      template <typename MapType>
      void load(const String& filename, MapType& map) throw (Exception::FileNotFound, Exception::ParseError)
      {
      	//try to open file
				QFile qfile(filename.c_str());
				if (!qfile.exists())
			    {
			      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
			    }

				QXmlSimpleReader parser;
				srand(static_cast<unsigned>(time(0)));
				parser.setFeature("http://xml.org/sax/features/namespaces",false);
				parser.setFeature("http://xml.org/sax/features/namespace-prefixes", false);

				map = MapType();  // clear map
				Internal::MzXMLHandler<MapType> handler(map);
				parser.setContentHandler(&handler);
				parser.setErrorHandler(&handler);

				QXmlInputSource source( &qfile );
				parser.parse(source);
      }

      /**
      	@brief Stores a map in a MzXML file.

      	@p map has to be a MSExperiment or have the same interface.
      */
      template <typename MapType>
      void store(const String& filename, const MapType& map)
			const throw (Exception::UnableToCreateFile)
			{
				std::ofstream os(filename.c_str());
				if (!os)
				{
					throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}

				//read data and close stream
				Internal::MzXMLHandler<MapType> handler(map);
				handler.writeTo(os);
				os.close();
			}
	};
} // namespace OpenMS

#endif // OPENMS_FOMAT_MZXMLFILE_H
