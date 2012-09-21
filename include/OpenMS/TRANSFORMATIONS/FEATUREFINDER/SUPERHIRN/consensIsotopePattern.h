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
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//

#ifndef CONSENS_ISOTOPE_PATTERN_H
#define CONSENS_ISOTOPE_PATTERN_H

#include <OpenMS/CONCEPT/Types.h>

#include <map>
#include <vector>

namespace OpenMS
{

class OPENMS_DLLAPI consensIsotopePattern{

  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  // stores teh consenus pattern:
  std::map<  double , double > isotopesTrace; 
  std::vector<  double > mzIsotopesStDev; 
  std::vector<  double > intensIsotopesStDev; 
  
  // stores the detected patterns by retention time
  std::map<double, std::pair< std::vector<double>, std::vector<double> > > rawIsotopes;
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
  
public:
    
    static double FT_MZ_TOLERANCE;
  
  // class destructor
  ~consensIsotopePattern();
  
  // class constructor
  consensIsotopePattern();
  // class copy constructor
  consensIsotopePattern(const consensIsotopePattern&);
  // class copy constructor
  consensIsotopePattern(const consensIsotopePattern*);
  
  
  //////////////////////////////////////////////////
  // overload operators:
  consensIsotopePattern& operator=(const consensIsotopePattern&);
  bool operator==(const consensIsotopePattern&);
  consensIsotopePattern& operator<=(const consensIsotopePattern&);
  consensIsotopePattern& operator>=(const consensIsotopePattern&);
  consensIsotopePattern& operator<(const consensIsotopePattern&);
  consensIsotopePattern& operator>(const consensIsotopePattern&);
  
  
  // construc the consenus pattern:
  void constructConsusPattern( );
  // order an isotope trace in the correct cluster:
  void addIsotopeTrace( double, double );
  // condens the pattern, make averge peaks from the traces:
  void condensIsotopePattern( std::pair< std::vector<double>, std::vector<double> >*);

  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class

  std::map<double, double>::iterator getConsensIsotopeIteratorStart(){return isotopesTrace.begin();}; 
  std::map<double, double>::iterator getConsensIsotopeIteratorEnd(){return isotopesTrace.end();}; 


};

} // ns

#endif

    
