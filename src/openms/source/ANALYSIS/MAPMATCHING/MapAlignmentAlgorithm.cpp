// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>

// Helper class to apply transformations to all sorts of KERNEL types
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>

// Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>

#include <OpenMS/KERNEL/ConversionHelper.h>

#include <OpenMS/CONCEPT/Factory.h>

using std::vector;

namespace OpenMS
{
  //register products here
  void MapAlignmentAlgorithm::registerChildren()
  {
    Factory<MapAlignmentAlgorithm>::registerProduct(MapAlignmentAlgorithmIdentification::getProductName(), &MapAlignmentAlgorithmIdentification::create);
    Factory<MapAlignmentAlgorithm>::registerProduct(MapAlignmentAlgorithmPoseClustering::getProductName(), &MapAlignmentAlgorithmPoseClustering::create);
    Factory<MapAlignmentAlgorithm>::registerProduct(MapAlignmentAlgorithmSpectrumAlignment::getProductName(), &MapAlignmentAlgorithmSpectrumAlignment::create);
  }

  MapAlignmentAlgorithm::MapAlignmentAlgorithm() :
    DefaultParamHandler("MapAlignmentAlgorithm"),
    ProgressLogger()
  {
  }

  MapAlignmentAlgorithm::~MapAlignmentAlgorithm()
  {
  }

  void MapAlignmentAlgorithm::alignPeakMaps(vector<MSExperiment<> > &, vector<TransformationDescription> &)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  void MapAlignmentAlgorithm::alignCompactFeatureMaps(vector<std::vector<Peak2D> > &, vector<TransformationDescription> &)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  void MapAlignmentAlgorithm::alignFeatureMaps(vector<FeatureMap> &, vector<TransformationDescription> &)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  void MapAlignmentAlgorithm::alignConsensusMaps(vector<ConsensusMap> & cms, vector<TransformationDescription> & tf)
  {
    LOG_WARN << "MapAlignmentAlgorithm::alignConsensusMaps() does not support ConsensusMaps directly. Converting to FeatureMaps.\n";

    vector<FeatureMap> maps_f;
    for (Size i = 0; i < cms.size(); ++i)
    {
      FeatureMap fm;
      MapConversion::convert(cms[i], true, fm);
      maps_f.push_back(fm);
    }
    // call FeatureMap version of group()
    alignFeatureMaps(maps_f, tf);
    // apply transform
    MapAlignmentTransformer::transformConsensusMaps(cms, tf);
  }

  void MapAlignmentAlgorithm::alignPeptideIdentifications(vector<vector<PeptideIdentification> > &, vector<TransformationDescription> &)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

  void MapAlignmentAlgorithm::setReference(Size reference_index,
                                           const String & reference_file)
  {
    if (reference_index || !reference_file.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "This algorithm does not support a reference for the alignment.");
    }
  }

  void MapAlignmentAlgorithm::fitModel(const String & model_type, const Param & params, vector<TransformationDescription> & trafos)
  {
    for (vector<TransformationDescription>::iterator it = trafos.begin();
         it != trafos.end(); ++it)
    {
      it->fitModel(model_type, params);
    }
  }

}
