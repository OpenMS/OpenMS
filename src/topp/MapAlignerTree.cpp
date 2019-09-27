// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Julia Thueringer $
// $Authors: Julia Thueringer $
// --------------------------------------------------------------------------

//#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTree.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_MapAlignerTree MapAlignerTree

  @brief Corrects retention time distortions between maps, using a tree and identifies features?
*/

  // We do not want this class to show up in the docu:
  /// @cond TOPPCLASSES

class TOPPMapAlignerTree :
  public TOPPMapAlignerBase
{

public:
  TOPPMapAlignerTree() :
    TOPPMapAlignerBase("MapAlignerTree", "Tree guided correction of retention time distortions between maps.")
  {}

private:
  template <typename MapType, typename FileType>
  void loadInputMaps_(vector<MapType>& maps, StringList& ins, FileType& input_file)
  {
      for (Size i = 0; i < ins.size(); ++i)
      {
        input_file.load(ins[i], maps[i]);
      }
  }

  template <typename MapType>
  void build_distance_matrix_(vector<MapType>& maps, vector<size_t>& dist_matrix)
  {
      for (Size i = 0; i < maps.size()-1; ++i)
      {
          for (Size j = i+1; j < maps.size(); ++j)
          {
              // get identified proteins of maps[i] and maps[j], sorted,
              // get union and intercept amount of proteins
              // create vectors for both maps containing RTs of identical proteins
              // pearsonCorrelationCoefficient(rt_map_i, rt_map_j)
              // dist_matrx[i][j] = pearson
          }
      }
  }

  void registerOptionsAndFlags_() override
  {
    String formats = "featureXML";
    TOPPMapAlignerBase::registerOptionsAndFlags_(formats, REF_NONE);
    //registerSubsection_("algorithm", "Algorithm parameters section");
    //registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
  }

  /*
  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "algorithm")
    {
      //MapAlignmentAlgorithmSpectrumAlignment algo;
      //return algo.getParameters();
    }
    if (section == "model")
    {
      //return getModelDefaults("interpolated");
    }
    //return Param(); // shouldn't happen
  }
  */

  ExitCodes main_(int, const char**) override
  {
    ExitCodes ret = checkParameters_();
    if (ret != EXECUTION_OK) return ret;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList in_files = getStringList_("in");
    StringList out_files = getStringList_("out");
    StringList out_trafos = getStringList_("trafo_out");


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    vector<FeatureMap> feature_maps(in_files.size());
    FeatureXMLFile fxml_file;
    loadMaps_(feature_maps, in_files, fxml_file);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    vector<size_t> dist_matrix;
    build_distance_matrix_(feature_maps, dist_matrix);

    // SingleLinkage tree = new SingleLinkage(dist_matrix)

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------


    return EXECUTION_OK;
  }

};

void same_feature_vecs_(std::vector< ProteinIdentification > vec1, std::vector< ProteinIdentification > vec2)
{

}

int main(int argc, const char** argv)
{
  TOPPMapAlignerTree tool;
  return tool.main(argc, argv);
}

/// @endcond
