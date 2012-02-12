// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:expandtab
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2011 -- Bastian Blank
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
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/CLUSTERING/HierarchicalClustering.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACPattern.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_SILACCLUSTERING_H
#define OPENMS_COMPARISON_CLUSTERING_SILACCLUSTERING_H

namespace OpenMS
{
  /**
   * @brief Clustering implementation for SILAC stuff.
   *
   * It cleans up the results of the hierarchical clustering for the search of
   * labeled and unlabeled peptide data. It removes too small clusters and joins
   * clusters with small gaps that the clustering split on purpose.
   *
   * @warning The cleanups will not work on random data
   */
  class OPENMS_DLLAPI SILACClustering : public HierarchicalClustering<SILACPattern *>
  {
    public:
      const DoubleReal rt_min;
      const DoubleReal rt_max_spacing;

      SILACClustering(const PointCoordinate &cluster_dimension, DoubleReal rt_min, DoubleReal rt_max_spacing)
        : HierarchicalClustering<SILACPattern *>(cluster_dimension), rt_min(rt_min), rt_max_spacing(rt_max_spacing)
      { }

      /**
       */
      void cluster();

    protected:
      /**
       * @brief Remove clusters smaller then rt_min
       */
      void removeSmall_();
      /**
       * @brief Join clusters with holes less then rt_max_spacing
       */
      void joinLarge_();
  };
}

#endif /* OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H */
