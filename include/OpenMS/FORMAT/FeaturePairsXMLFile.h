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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FEATUREPAIRSXMLFILE_H
#define OPENMS_FORMAT_FEATUREPAIRSXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/FeaturePairsHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

  /**
    @brief This class provides Input/Output functionality for the class DFeaturePairVector.
    
    The features pairs are computed by an instance of DBaseFeatureMatcher during the
    matching of MS maps. The features pairs are stored in a pseudo XML format. No schema
    has been developed yet therefore no validation can be performed.

    @ingroup FileIO
  */
  class FeaturePairsXMLFile
        : public Internal::XMLFile
  {
	  public:
	    /** @name Constructors and Destructor */
	    //@{
	    ///Default constructor
	    FeaturePairsXMLFile();
	    ///Destructor
	    ~FeaturePairsXMLFile();
	    //@}

	    /// loads the file with name @p filename into @p pairs.
	    void load(String filename, std::vector< ElementPair < Feature > > & pairs) throw (Exception::FileNotFound, Exception::ParseError);
	
	    /// stores the pair vector @p pairs in file with name @p filename.
	    void store(String filename, const std::vector< ElementPair < Feature > > & pairs) const throw (Exception::UnableToCreateFile);
	
	    /**
	      @brief Convert pair vector into feature map
	
	    */
	    static void pairsToFeatures(const std::vector< ElementPair < Feature > >& pairs, FeatureMap<>& map);
	
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_FEATUREPAIRSXMLFILE_H
