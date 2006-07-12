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

#include <qxml.h>

#include <iostream>  
#include <fstream>   
#include <sstream>
#include <stdio.h>
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
      DFeaturePairsFile() { initParser_(); }
       ///Destructor
      virtual ~DFeaturePairsFile() { delete parser_; }
      //@}

      /** @name Accessors */
      //@{
      /// loads the file with name @p filename into @p pairs.
      template<Size D> 
      void load(String filename, DFeaturePairVector<D>& pairs) throw (Exception::FileNotFound)
      {
				Internal::DFeaturePairsHandler<D> handler(pairs);
				parser_->setContentHandler(&handler);
				parser_->setErrorHandler(&handler);

				QFile qfile(filename.c_str());
				if (!qfile.exists())
				{	
				  throw new Exception::FileNotFound(__FILE__,__LINE__,"DFeaturePairsFile::load()",filename);
				} 	
				QXmlInputSource source(qfile);
				parser_->parse(source);
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

				Internal::DFeaturePairsHandler<D> handler(pairs);
				handler.writeTo(os);
				os.close();
      }
               
      //@}
      
     protected:
      /// parser initialization
			inline void initParser_()
			{
				srand(static_cast<unsigned>(time(0)));
				parser_ = new QXmlSimpleReader();
				parser_->setFeature("http://xml.org/sax/features/namespaces",false);
				parser_->setFeature("http://xml.org/sax/features/namespace-prefixes", false);
			}
			
			/// the parser
			QXmlSimpleReader* parser_;
			
	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_DFEATUREPAIRSFILE_H
