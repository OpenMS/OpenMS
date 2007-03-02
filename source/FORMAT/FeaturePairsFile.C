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

#include<OpenMS/FORMAT/FeaturePairsFile.h>

namespace OpenMS
{

  // not much going on here...

  FeaturePairsFile::FeaturePairsFile()
      : SchemaFile()
  {
  }
  ///Destructor
  FeaturePairsFile::~FeaturePairsFile()
  {
  }

  void FeaturePairsFile::pairsToFeatures(const std::vector< ElementPair < Feature > > & pairs, FeatureMap<>& map)
  {
    map.clear();

    for (Size i=0; i<pairs.size(); ++i)
    {
      map.push_back(pairs[i].getFirst());
      map.push_back(pairs[i].getSecond());
    }
  }

  void FeaturePairsFile::load(String filename, std::vector< ElementPair < Feature > > & pairs) throw (Exception::FileNotFound, Exception::ParseError)
  {
    Internal::FeaturePairsHandler handler(pairs,filename);
    parse_(filename, &handler);
  }

  /// stores the pair vector @p pairs in file with name @p filename.
  void FeaturePairsFile::store(String filename, const std::vector< ElementPair < Feature > > & pairs) const throw (Exception::UnableToCreateFile)
  {
    if (pairs.empty()) return;
    Internal::FeaturePairsHandler handler(pairs,filename);
    save_(filename, &handler);
  }

}
