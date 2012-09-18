// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
