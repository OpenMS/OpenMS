// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//
//  by Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#ifndef _USE_EXTERNAL_ISOTOPIC_DISTRIBUTION_H
#define _USE_EXTERNAL_ISOTOPIC_DISTRIBUTION_H

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{

class OPENMS_DLLAPI ExternalIsotopicDistribution{

  
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:

  static std::multimap< double, PeptideIsotopeDisribution> allExternalPepdistributions;
  
  ////////////////////////////////////////////////
  // declaration of the public members:

public:
  
  static bool EXTERNAL_ISOTOPIC_PROFILES;  
  static std::string XMLInputFile;
  static double EXTERNAL_DISTRIBUTION_MONO_ISOTOPE_PPM_TOLERANCE;

  // class destructor
  ~ExternalIsotopicDistribution();
  
  // class constructor
  ExternalIsotopicDistribution();
  // class copy constructor
  ExternalIsotopicDistribution(const ExternalIsotopicDistribution&);
  // class copy constructor
  ExternalIsotopicDistribution(const ExternalIsotopicDistribution*);
  
  
  //////////////////////////////////////////////////
  // overload operators:
  ExternalIsotopicDistribution& operator=(const ExternalIsotopicDistribution&);
  bool operator==(const ExternalIsotopicDistribution&);
  ExternalIsotopicDistribution& operator<=(const ExternalIsotopicDistribution&);
  ExternalIsotopicDistribution& operator>=(const ExternalIsotopicDistribution&);
  ExternalIsotopicDistribution& operator<(const ExternalIsotopicDistribution&);
  ExternalIsotopicDistribution& operator>(const ExternalIsotopicDistribution&);
  
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  
  
  
  
  //////////////////////////////////////////////////
  // function to extract external isotopic profiles
  // by an input monoisotopic mass

  static PeptideIsotopeDisribution* extractExternalIsotopicProfile( double, int , double);
  //  check if two masses to be the same
  static bool checkMonoIsotopOverlap( double , int, double, PeptideIsotopeDisribution* );
  
  //////////////////////////////////////////////////
  // parse the external isotopic profiles
  // from xml file
  //  static void parseExternalIsotopicProfiles( );
  // extract isotopic profiles
  //  static void extractIsotopicProfiles( TiXmlDocument* );

  // init the retention time segments
  static void initRetentionTimeSegments( double, double);

  
};

} // ns

#endif
