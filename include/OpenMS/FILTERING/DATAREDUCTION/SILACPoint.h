// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut, Steffen Sass, Bastian Blank $
// --------------------------------------------------------------------------


#ifndef OPENMS_FILTERING_DATAREDUCTION_SILACPOINT_H
#define OPENMS_FILTERING_DATAREDUCTION_SILACPOINT_H

namespace OpenMS
{
  /**
   * @brief A single SILAC point
   * @see HashGrid
   * @ingroup Datastructures
   */
  class OPENMS_DLLAPI SILACPoint
  {
    public:
      /**
       * @brief m/z value of the element
       */
      DoubleReal mz;

      /**
       *@brief RT value of the element
       */
      DoubleReal rt;

      /**
       * @brief exact m/z positions at which the intensities are read
       */
      std::vector<std::vector<DoubleReal> > mz_positions;
      
      /**
       * @brief intensities at RT and the exact m/z positions
       */
      std::vector<std::vector<DoubleReal> > intensities;
      
      /**
       * @brief mass shifts [Da] used in the filter
       */	
      std::vector<DoubleReal> mass_shifts;

      /**
       * @brief charge of the cluster (i.e. peptide) which the data point is part of
       */
      Int charge;

      /**
       * @brief number of isotopes per peptide of the cluster
       */
      Int isotopes_per_peptide;

      /**
       * @brief quality of the cluster
       */
      DoubleReal quality;

      SILACPoint()
        : mz(0), rt(0), charge(0), isotopes_per_peptide(0), quality(0)
      { }
  };
}

#endif
