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
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

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

  /// Type to store retention times given for individual peptide sequence
  typedef std::map<String, double> SeqAndRT;

  template <typename MapType, typename FileType>
  void loadInputMaps_(vector<MapType>& maps, StringList& ins, FileType& input_file)
  {
      for (Size i = 0; i < ins.size(); ++i)
      {
        input_file.load(ins[i], maps[i]);
      }
  }

  void getPeptideSequences_(vector<PeptideIdentification>& peptides, SeqAndRT& peptide_rts)
  {
      for (std::vector<PeptideIdentification>::iterator pep_it =
        peptides.begin(); pep_it != peptides.end(); ++pep_it)
      {
        if (!pep_it->getHits().empty())
        {
          const String& sequence = pep_it->getHits()[0].getSequence().toString();
          peptide_rts[sequence] = pep_it->getRT();
        }
      }
  }

  template <typename MapType>
  void build_distance_matrix_(Size maps_amount, vector<SeqAndRT>& maps_peptides, vector<size_t>& dist_matrix)
  {
      /*
      for (Size i = 0; i < maps_amount.size()-1; ++i)
      {
          for (Size j = i+1; j < maps_amount.size(); ++j)
          {
              // get identified proteins of maps[i] and maps[j], sorted,
              // getRetentionTimes_(maps[i], maps_rts[i]);
              // get union and intercept amount of proteins
              // create vectors for both maps containing RTs of identical proteins
              // pearsonCorrelationCoefficient(rt_map_i, rt_map_j)
              // dist_matrx[i][j] = pearson
          }
      }
      */
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
    loadInputMaps_(feature_maps, in_files, fxml_file);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    // one set of RT data for each input map
    vector<SeqAndRT> maps_peptides(in_files.size());
    // get Peptide/ RT tuple for all features, seperated by input file
    for (vector<FeatureMap>::iterator maps_it = feature_maps.begin(); maps_it != feature_maps.end(); ++maps_it)
    {
      for (vector<Feature>::iterator feature_it = maps_it->begin();
         feature_it != maps_it->end(); ++feature_it)
      {
        if (feature_it->getPeptideIdentifications().size()>0)
        {
            getPeptideSequences_(feature_it->getPeptideIdentifications(), maps_peptides[distance(feature_maps.begin(), maps_it)]);
        }
      }
    }
    /* for debug
    for (vector<SeqAndRT>::iterator itmaps= maps_peptides.begin(); itmaps != maps_peptides.end(); ++itmaps)
    {
        for (SeqAndRT::iterator it = itmaps->begin(); it != itmaps->end(); ++it)
        {
            std::cout << it->first << " : " << it->second << "\n";
        }
    }*/

    // remove feature-maps? or is there another way to get relevant data for later alignment?


    // RTs of cluster (only petides present in both parent maps)
    //vector<SeqAndRT> clusters_rts(in_files.size()*in_files.size());
    vector<size_t> dist_matrix(in_files.size()*in_files.size());
    build_distance_matrix_(feature_maps.size(), maps_peptides, dist_matrix);

    // SingleLinkage tree = new SingleLinkage(dist_matrix)

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------


    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPMapAlignerTree tool;
  return tool.main(argc, argv);
}

/// @endcond
