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

#ifndef OPENMS_FORMAT_DGRIDFILE_H
#define OPENMS_FORMAT_DGRIDFILE_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/DGridHandler.h>

#include <qxml.h>

#include <iostream.h>  
#include <fstream.h>   
#include <sstream>
#include <stdio.h>
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
      DGridFile() { initParser_(); }
       ///Destructor
      virtual ~DGridFile() { delete parser_; }
      //@}

      /** @name Accessors */
      //@{
      /// loads the file with name @p filename into @p grid. 
      template<Size D> 
      void load(String filename, DGrid<D>& grid) throw (Exception::FileNotFound)
      {
				Internal::DGridHandler<D> handler(grid);
				parser_->setContentHandler(&handler);
				parser_->setErrorHandler(&handler);

				QFile qfile(filename.c_str());
				if (!qfile.exists())
				{	
				  throw new Exception::FileNotFound(__FILE__,__LINE__,"DGridFile::load()",filename);
				} 	
				QXmlInputSource source(qfile);
				parser_->parse(source);
      }

      /// stores the grid @p grid in file with name @p filename. 
      template<Size D> 
      void store(String filename, const DGrid<D>& grid) const throw (Exception::UnableToCreateFile)
      {
			  if (grid.empty()) return;

				ofstream os(filename.c_str(),fstream::out);
				if (!os.is_open())
				{
					os.close();
					throw Exception::UnableToCreateFile(__FILE__,__LINE__,"DGridFile::store()",filename);
				}

				Internal::DGridHandler<D> handler(grid);
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

#endif // OPENMS_FORMAT_DGRIDFILE_H
