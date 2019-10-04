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
#include <OpenMS/KERNEL/ConversionHelper.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>


using namespace OpenMS;
using namespace Math;
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

  class PeptideIdentificationsPearsonDistance
  {
  public:
    /*
    // default constructor
    PeptideIdentificationsPearsonDistance();

    // copy constructor
    PeptideIdentificationsPearsonDistance(const PeptideIdentificationsPearsonDistance & source);

    // destructor
    ~PeptideIdentificationsPearsonDistance();
    // assignment operator
    PeptideIdentificationsPearsonDistance & operator=(const PeptideIdentificationsPearsonDistance & source)
    {
      if (this != &source)
      {
        DefaultParamHandler::operator=(source);
      }
      return *this;
    }
    */

    float operator()(const SeqAndRT& map_first, const SeqAndRT& map_second) const
    {
      // create vectors for both maps containing RTs of identical proteins and
      // get union and intercept amount of proteins
      SeqAndRT::const_iterator pep1_it = map_first.begin();
      SeqAndRT::const_iterator pep2_it = map_second.begin();
      vector<double> intercept_rts1;
      vector<double> intercept_rts2;
      while (pep1_it != map_first.end() && pep2_it != map_second.end())
      {
          if (pep1_it->first < pep2_it->first)
          {
              ++pep1_it;
          }
          else if (pep2_it->first < pep1_it->first)
          {
              ++pep2_it;
          }
          else
          {
              //intercept_peps1[pep1_it->first] = pep1_it->second; //if second DoubleList, then pushback of both possible
              //intercept_peps2[pep2_it->first] = pep2_it->second;
              intercept_rts1.push_back(pep1_it->second);
              intercept_rts2.push_back(pep2_it->second);
              ++pep1_it;
              ++pep2_it;
          }
      }
      Size intercept_size = intercept_rts1.size();
      SeqAndRT union_map_tmp;
      union_map_tmp.insert(map_first.begin(), map_first.end());
      union_map_tmp.insert(map_second.begin(), map_second.end());
      Size union_size = union_map_tmp.size();

      // pearsonCorrelationCoefficient(rt_map_i, rt_map_j)
      float pearson_val = static_cast<float>(pearsonCorrelationCoefficient(intercept_rts1.begin(), intercept_rts1.end(), intercept_rts2.begin(), intercept_rts2.end()));

      if (pearson_val > 1)
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
      return (1 - (pearson_val*intercept_size/union_size));
    }

    /*
    // calculates self similarity
    virtual float operator()(const PeptideIdentificationsPearsonDistance & a) const = 0;

    // registers all derived products
    static void registerChildren();

    static const String getProductName()
    {
        return "PeptideIdentificationsPearsonDistance";
    }
    */

  }; // end of PeptidIdentificationsPearsonDifferencer

  template <typename MapType, typename FileType>
  void loadInputMaps_(vector<MapType>& maps, StringList& ins, FileType& input_file)
  {
      ProgressLogger progresslogger;
      progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
      progresslogger.startProgress(0, ins.size(), "loading input files");
      for (Size i = 0; i < ins.size(); ++i)
      {
        progresslogger.setProgress(i);
        input_file.load(ins[i], maps[i]);
      }
      progresslogger.endProgress();
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

  /*
  void build_distance_matrix_(Size maps_amount, vector<SeqAndRT>& maps_peptides, DistanceMatrix<float>& dist_matrix)
  {
      for (Size i = 0; i < maps_amount-1; ++i)
      {
          for (Size j = i+1; j < maps_amount; ++j)
          {
              // create vectors for both maps containing RTs of identical proteins and
              // get union and intercept amount of proteins
              SeqAndRT::iterator pep1_it = maps_peptides[i].begin();
              SeqAndRT::iterator pep2_it = maps_peptides[j].begin();
              vector<double> intercept_rts1;
              vector<double> intercept_rts2;
              while (pep1_it != maps_peptides[i].end() && pep2_it != maps_peptides[j].end())
              {
                  if (pep1_it->first < pep2_it->first)
                  {
                      ++pep1_it;
                  }
                  else if (pep2_it->first < pep1_it->first)
                  {
                      ++pep2_it;
                  }
                  else
                  {
                      //intercept_peps1[pep1_it->first] = pep1_it->second; //if second DoubleList, then pushback of both possible
                      //intercept_peps2[pep2_it->first] = pep2_it->second;
                      intercept_rts1.push_back(pep1_it->second);
                      intercept_rts2.push_back(pep2_it->second);
                      ++pep1_it;
                      ++pep2_it;
                  }
              }
              Size intercept_size = intercept_rts1.size();
              SeqAndRT union_map_tmp;
              union_map_tmp.insert(maps_peptides[i].begin(), maps_peptides[i].end());
              union_map_tmp.insert(maps_peptides[j].begin(), maps_peptides[j].end());
              Size union_size = union_map_tmp.size();

              // pearsonCorrelationCoefficient(rt_map_i, rt_map_j)
              float pearson_val = static_cast<float>(pearsonCorrelationCoefficient(intercept_rts1.begin(), intercept_rts1.end(), intercept_rts2.begin(), intercept_rts2.end()));

              dist_matrix.setValue(i,j, 1-pearson_val*intercept_size/union_size);
          }
      }
  }
  */

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
    Size in_files_size = in_files.size();
    vector<FeatureMap> feature_maps(in_files_size);
    FeatureXMLFile fxml_file;
    loadInputMaps_(feature_maps, in_files, fxml_file);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    // one set of RT data for each input map
    vector<SeqAndRT> maps_peptides(in_files_size);
    // get Peptide/ RT tuple for all features, seperated by input file
    for (vector<FeatureMap>::iterator maps_it = feature_maps.begin(); maps_it != feature_maps.end(); ++maps_it)
    {
      for (vector<Feature>::iterator feature_it = maps_it->begin();
         feature_it != maps_it->end(); ++feature_it)
      {
        if (feature_it->getPeptideIdentifications().size()>0)
        {
            getPeptideSequences_(feature_it->getPeptideIdentifications(), maps_peptides[static_cast<unsigned long>(distance(feature_maps.begin(), maps_it))]);
        }
      }
    }

    // remove feature-maps? or is there another way to get relevant data for later alignment?


    // RTs of cluster (only petides present in both parent maps)
    //vector<SeqAndRT> clusters_rts(in_files.size()*in_files.size());
    //DistanceMatrix<float> dist_matrix(in_files_size*(in_files_size-1)/2);
    //build_distance_matrix_(feature_maps.size(), maps_peptides, dist_matrix);

    PeptideIdentificationsPearsonDistance pepDist;
    SingleLinkage sl;
    std::vector<BinaryTreeNode> tree;
    DistanceMatrix<float> dist_matrix;
    ClusterHierarchical ch;
    ch.cluster<SeqAndRT, PeptideIdentificationsPearsonDistance>(maps_peptides, pepDist, sl, tree, dist_matrix);

    // to print tree
    //ClusterAnalyzer ca;
    //std::cout << ca.newickTree(tree) << std::endl;



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
