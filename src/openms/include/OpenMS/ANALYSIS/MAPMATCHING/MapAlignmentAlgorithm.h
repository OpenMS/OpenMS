// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHM_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHM_H

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  /**
      @brief Base class for all map-alignment algorithms

      It takes two or more maps and corrects for retention time distortions.

      The input maps are transformed and the transformation description is returned.

      @improvement The maps should not be all loaded before the algorithm  - in order to save memory e.g. in the star-wise approach (Clemens)
  */
  class OPENMS_DLLAPI MapAlignmentAlgorithm :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    MapAlignmentAlgorithm();

    /// Destructor
    virtual ~MapAlignmentAlgorithm();

    /**
    @brief Aligns peak maps

    @exception Exception::NotImplemented is thrown if an algorithm cannot align peak maps
    */
    virtual void alignPeakMaps(std::vector<MSExperiment<> >&, std::vector<TransformationDescription>&);

    /**
    @brief Aligns vectors of 2D peaks (memory efficient version of FeatureMap)

    @exception Exception::NotImplemented is thrown if an algorithm cannot align feature maps
    */
    virtual void alignCompactFeatureMaps(std::vector<std::vector<Peak2D> >&, std::vector<TransformationDescription>&);

    /**
    @brief Aligns feature maps

    @exception Exception::NotImplemented is thrown if an algorithm cannot align feature maps
    */
    virtual void alignFeatureMaps(std::vector<FeatureMap>&, std::vector<TransformationDescription>&);

    /**
    @brief Aligns consensus maps

    @exception Exception::NotImplemented is thrown if an algorithm cannot align consensus maps
    */
    virtual void alignConsensusMaps(std::vector<ConsensusMap>&, std::vector<TransformationDescription>&);

    /**
    @brief Aligns peptide identifications

    @exception Exception::NotImplemented is thrown if an algorithm cannot align peptide identifications
    */
    virtual void alignPeptideIdentifications(std::vector<std::vector<PeptideIdentification> >&, std::vector<TransformationDescription>&);

    /**
         @brief Defines a reference for the alignment

         @exception Exception::NotImplemented The algorithm does not support references
        */
    template <typename MapType>
    virtual void setReference(const MapType& map)
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

    /**
         @brief Fits a model with given parameters to the transformations

         This will not alter transformations of reference files (transformation type "identity").
    */
    static void fitModel(const String& model_type, const Param& params, std::vector<TransformationDescription>& trafos);

    /// Register all derived classes in this method
    static void registerChildren();

private:
    /// Copy constructor is not implemented -> private
    MapAlignmentAlgorithm(const MapAlignmentAlgorithm&);
    /// Assignment operator is not implemented -> private
    MapAlignmentAlgorithm& operator=(const MapAlignmentAlgorithm&);

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHM_H
