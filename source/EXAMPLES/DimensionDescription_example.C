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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------


#include <OpenMS/KERNEL/DimensionDescription.h>

#include <iostream>


template < unsigned int D > struct SomeClass { unsigned int getD () { return D; } };
template < unsigned int D > unsigned int getD () { return D; }

int
main()
{
  
  using namespace OpenMS;

  using std::cout;
  using std::endl;



  // Explaining the general technique
  cout << "Using 17 as template arg" << endl;
  cout << "from class:     " << SomeClass<17>().getD() << endl;
  cout << "from function:  " << getD<17>() << endl;



  // Now we exemplify the DimensionDescription class template
  typedef DimensionDescription < LCMS_Tag > DimDesc;
  typedef DimDesc::DimensionId DimensionId;  


  DimensionId dimension_id = DimDesc::MZ;
  cout << "MZ: " << dimension_id << endl;

  dimension_id = DimDesc::RT;
  cout << "RT: " << dimension_id << endl;
  
  cout << "Using MZ as template arg" << endl;
  cout << "from class:     " << SomeClass<DimDesc::MZ>().getD() << endl;
  cout << "from function:  " << getD<DimDesc::MZ>() << endl;

  cout << "Using RT as template arg" << endl;
  cout << "from class:     " << SomeClass<DimDesc::RT>().getD() << endl;
  cout << "from function:  " << getD<DimDesc::RT>() << endl;


  // wow! even this compiles:
  cout << "Using const_dimension_id == RT as template arg" << endl;
  const DimensionId const_dimension_id = DimDesc::RT;
  cout << "from class:     " << SomeClass<const_dimension_id>().getD() << endl;
  cout << "from function:  " << getD<const_dimension_id>() << endl;


  cout << "Now here is the information about each dimension...\n";
  for ( int dim = 0; dim < DimDesc::DIMENSION; ++dim )
	{
		cout
			<< dim << ' '
			<< DimDesc::dimension_name_short[dim] << ' '
			<< DimDesc::dimension_name_full[dim] << ' '
			<< DimDesc::dimension_unit_short[dim] << ' '
			<< DimDesc::dimension_unit_full[dim] << ' '
			<< endl;
	}


  // Of course, you don't need to prefix "DimensionDescription < LCMS_Tag >::" all the time
  unsigned int const & MZ = DimDesc::MZ;
  unsigned int const & RT = DimDesc::RT;

  cout
    << "MZ="   << DimDesc::dimension_name_short[MZ]
    << "\nRT=" << DimDesc::dimension_name_short[RT]
    << endl;

  return 0;
}
