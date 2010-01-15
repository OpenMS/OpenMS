// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGAFFINESUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGAFFINESUPERIMPOSER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>

namespace OpenMS
{

  /**
   @brief A superimposer that uses a voting scheme, also known as pose clustering,
   to find a good affine transformation.

   This algorithm works on two consensus maps.  It computes an affine
   transformation that maps the elements of second map as near as possible
   to the elements in the first map.

   The voting scheme hashes affine transformations between pairs of features in
   map one and pairs of features in map two.  Each such pair of pairs defines a
   (potential) "pose" of the second map relative to the first.
   Then it finds a cluster in the parameter space of these poses.
   The affine transformation is then computed from this
   cluster of potential poses, hence the name pose clustering.

   @sa PoseClusteringShiftSuperimposer

   @htmlinclude OpenMS_PoseClusteringAffineSuperimposer.parameters

   @ingroup MapAlignment
   */
  class OPENMS_DLLAPI PoseClusteringAffineSuperimposer : public BaseSuperimposer
  {
  public:

    /// Default ctor
    PoseClusteringAffineSuperimposer();

    /// Destructor
    virtual
    ~PoseClusteringAffineSuperimposer()
    {
    }

    /**
     @brief Estimates the transformation and fills the given mapping function. (Has a precondition!)

     @note Exactly two input maps must be given.

     @pre  For performance reasons, we trust that (the equivalent of:) <code>maps[0].updateRanges(); maps[1].updateRanges();</code> has been done <i>before</i> calling this.  You have been warned!

     @exception IllegalArgument is thrown if the input maps are invalid.
     */
    virtual void
    run(const std::vector<ConsensusMap>& maps, std::vector<TransformationDescription>& transformations);

    /// Returns an instance of this class
    static BaseSuperimposer*
    create()
    {
      return new PoseClusteringAffineSuperimposer();
    }

    /// Returns the name of this module
    static const String
    getProductName()
    {
      return "poseclustering_affine";
    }

  };

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGAFFINESUPERIMPOSER_H
