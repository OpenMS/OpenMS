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

#ifndef OPENMS_FORMAT_DFEATUREMAPFILE_H
#define OPENMS_FORMAT_DFEATUREMAPFILE_H

#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/DFeatureMapHandler.h>

#include <qxml.h>

#include <iostream.h>
#include <fstream.h>
#include <sstream>
#include <stdio.h>

namespace OpenMS
{
	
	/**
  	@brief This class provides Input/Output functionality for DFeatureMaps
  		
  		
  
  	@ingroup FileIO
  */
  class DFeatureMapFile
  {
	 public:
		/** @name Constructors and Destructor */
		//@{

		///Default constructor
		DFeatureMapFile() { initParser_(); }
		///Destructor
		virtual ~DFeatureMapFile() { delete parser_; }
		//@}

		/** @name Accessors */
		//@{
		/// loads the file with name @p filename into @p map. General case is not implemented!
		template<Size D> 
		void load(String filename, DFeatureMap<D>& map) throw (Exception::FileNotFound);
			
		
		/// loads the file with name @p filename into @p map.
		void load(String filename, DFeatureMap<2>& feature_map) throw (Exception::FileNotFound)
		{
			feature_map.clear();

			Internal::DFeatureMapHandler<2> handler(feature_map);
			parser_->setContentHandler(&handler);
			parser_->setErrorHandler(&handler);
		
			QFile qfile(filename.c_str());
			if (!qfile.exists())
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, "DFeatureMapFile::load",filename);
			}
			QXmlInputSource source( qfile );
			parser_->parse(source);
		}
      
		/// stores the map @p map in file with name @p filename. General case is not implemented!
		template<Size D> 
		void store(String filename, const DFeatureMap<D>& map) const throw (Exception::UnableToCreateFile);
			
		
		/// stores the map @p feature_map in file with name @p filename.
		void store(String filename, const DFeatureMap<2>& feature_map) const throw (Exception::UnableToCreateFile)
		{
			if (feature_map.empty()) return;
      
			// open file
			ofstream os(filename.c_str(), fstream::out);
			if (!os.is_open())
			{
				os.close();
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, "DFeatureMapFile::store",filename);
			}	
    		
			//read data and close stream
			Internal::DFeatureMapHandler<2> handler(feature_map);
			handler.writeTo(os);
			os.close();
		}
                
		//@}
      
	 protected:
		// parser initialization
		inline void initParser_()
		{
			srand(static_cast<unsigned>(time(0)));
			parser_ = new QXmlSimpleReader();
			parser_->setFeature("http://xml.org/sax/features/namespaces",false);
			parser_->setFeature("http://xml.org/sax/features/namespace-prefixes", false);
		}
		// the parser
		QXmlSimpleReader* parser_;

	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_DEFEATUREMAPFILE_H
