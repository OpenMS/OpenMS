// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>

#include "../topp/FeatureLinkerBase.cpp"

#include <fstream>
#include <iostream>

using namespace OpenMS;
using namespace std;

//
//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_FeatureLinkerUnlabeledKD FeatureLinkerUnlabeledKD

 @brief Group corresponding features across labelfree experiments.

 Group corresponding features across labelfree experiments. This tool
 produces results similar to those of FeatureLinkerUnlabeledQT, since it
 optimizes a similar objective. However, this algorithm is more efficient
 than FLQT as it uses a kd-tree for fast 2D region queries in m/z - RT space
 and a sorted binary search tree to choose the best cluster among the remaining
 ones in O(1). Insertion and searching in this tree have O(log n) runtime.
 KD-tree insertion and search have O(log n) runtime. The overall complexity of
 the algorithm is O(n log(n)) time and O(n) space.

 In practice, the runtime of FeatureLinkerUnlabeledQT is often not
 significantly worse than that of FeatureLinkerUnlabeledKD if the datasets
 are relatively small and/or the value of the -nr_partitions parameter is
 chosen large enough. If, however, the datasets are very large, and especially
 if they are so dense that a partitioning based on the specified m/z
 tolerance is not possible anymore, then this algorithm becomes orders of
 magnitudes faster than FLQT.

 Notably, this algorithm can be used to align featureXML files containing
 unassembled mass traces (as produced by MassTraceExtractor), which is often
 impossible for reasonably large datasets using other aligners, as these
 datasets tend to be too dense and hence cannot be partitioned.

 Prior to feature linking, this tool performs an (optional) retention time
 transformation on the features using LOWESS regression in order to minimize
 retention time differences between corresponding features across different
 maps. These transformed RTs are used only internally. In the results, original
 RTs will be reported.

 The linking behavior can be influenced by separately specifying how to use
 the available charge and adduct information.  Options  allow to restrict
 linking to features with the same adduct/charge (or lack thereof, i.e.
 features with charge zero or no adduct annotation), additionally allowing
 the linking of charged/adduct-annotated features with those having no
 charge/adduct information, or allowing all features to be linked irrespective
 of charge state/adduct information.
 
 Note that the more relaxed the allowed grouping criteria, the larger internally
 used connected components memory-wise. More stringent m/z or retention time
 tolerances might be required then.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_FeatureLinkerUnlabeledKD.cli
 <B>INI file documentation of this tool:</B>
 @htmlinclude TOPP_FeatureLinkerUnlabeledKD.html

 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureLinkerUnlabeledKD :
  public TOPPFeatureLinkerBase
{

public:
  TOPPFeatureLinkerUnlabeledKD() :
    TOPPFeatureLinkerBase("FeatureLinkerUnlabeledKD", "Groups corresponding features from multiple maps.")
  {
    setLogType(CMD);
  }

protected:
  void registerOptionsAndFlags_() override
  {
    TOPPFeatureLinkerBase::registerOptionsAndFlags_();
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    FeatureGroupingAlgorithmKD algo;
    Param p = algo.getParameters();
    return p;
  }

  ExitCodes main_(int, const char **) override
  {
    FeatureGroupingAlgorithmKD algo;
    return TOPPFeatureLinkerBase::common_main_(&algo);
  }

};


int main(int argc, const char ** argv)
{
  TOPPFeatureLinkerUnlabeledKD tool;
  return tool.main(argc, argv);
}

/// @endcond
