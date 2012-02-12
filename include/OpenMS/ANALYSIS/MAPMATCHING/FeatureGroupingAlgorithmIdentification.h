// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMIDENTIFICATION_H
#define OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMIDENTIFICATION_H

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

namespace OpenMS
{
  /**
	 @deprecated Deprecated in OpenMS 1.7.

   @brief A map feature grouping algorithm for identified features.

   It takes many maps and searches for corresponding features.
   The corresponding features must be aligned, but may have small position deviations.

   @htmlinclude OpenMS_FeatureGroupingAlgorithmIdentification.parameters

   @ingroup FeatureGrouping
   */
  class OPENMS_DLLAPI FeatureGroupingAlgorithmIdentification : public FeatureGroupingAlgorithm
  {
  public:
    /// Default constructor
    FeatureGroupingAlgorithmIdentification();

    /// Destructor
    virtual
    ~FeatureGroupingAlgorithmIdentification();

    /**
     @brief Applies the algorithm

     @exception IllegalArgument is thrown if less than two input maps are given.
     */
    virtual void
    group(const std::vector<FeatureMap<> >& maps, ConsensusMap& out);

    /// Creates a new instance of this class (for Factory)
    static FeatureGroupingAlgorithm*
    create()
    {
      return new FeatureGroupingAlgorithmIdentification();
    }

    /// Returns the product name (for the Factory)
    static String
    getProductName()
    {
      return "identification";
    }

  private:

    /// Copy constructor intentionally not implemented -> private
    FeatureGroupingAlgorithmIdentification(const FeatureGroupingAlgorithmIdentification&);
    /// Assignment operator intentionally not implemented -> private
    FeatureGroupingAlgorithmIdentification&
    operator=(const FeatureGroupingAlgorithmIdentification&);

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMIDENTIFICATION_H
