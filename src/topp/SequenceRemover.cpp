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
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;
using namespace boost;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_SequenceRemover SequenceRemover

    @brief Removes 'n' sequences from peptides randomly from featureXML files.

  Input and output format are 'featureXML'. The tool allows you to remove known sequences
  from peptides in featureXMLs.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSequenceRemover :
        public TOPPMapAlignerBase
{
public:
  TOPPSequenceRemover() :
          TOPPMapAlignerBase("SequenceRemover", "Removes features from featureXML files if they contain randomly selected sequences.")
  {

  }

  // function declarations
  void loadInputMap_(FeatureMap& map, const String& in, FeatureXMLFile& fxml_file);
  void storeOutFiles_(const vector<FeatureMap>& out_maps, const StringList& out_files, FeatureXMLFile& fxml_file);
  void storeRemovedOutFiles_(const vector<FeatureMap>& removed_maps, const StringList& removed_out_files, FeatureXMLFile& fxml_file);

protected:
  void registerOptionsAndFlags_() override
  {
    TOPPMapAlignerBase::registerOptionsAndFlags_("featureXML", REF_NONE);
    registerOutputFileList_("removed_out", "<files>", StringList(), "Output files containing removed features.", true);
    setValidFormats_("removed_out", {"featureXML"});
    registerDoubleOption_("percent_to_remove", "<double>", 0.1, "Percentage of peptide identifications to be remove", false, false);
    setMinFloat_("percent_to_remove", 0.0);
    setMaxFloat_("percent_to_remove", 1.0);
    registerIntOption_("map_index", "<int>", 1, "Index of featureXML file from which peptide identifications are randomly chosen; starts with 1", false, false);
    setMinInt_("map_index", 1);
  }

  ExitCodes main_(int, const char**) override
  {
    FeatureXMLFile fxml_file;
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList ins = getStringList_("in");
    StringList outs = getStringList_("out");
    StringList removed_outs = getStringList_("removed_out");
    double percent_to_remove = getDoubleOption_("percent_to_remove");
    Size map_index = getIntOption_("map_index")-1;

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    vector<FeatureMap> maps(ins.size());

    if (map_index > (ins.size()-1))
    {
      OPENMS_LOG_WARN << "Map_index doesn't match the number of input files. Choose peptide Identifications from input file 1" << endl;
      map_index = 0;
    }

    loadInputMap_(maps[map_index], ins[map_index], fxml_file);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    // determine amount of peptide ids to be removed
    Size to_be_removed = (Size)ceil(maps[map_index].size()*percent_to_remove);
    OPENMS_LOG_INFO << to_be_removed << "peptide identifications to be removed." << std::endl;

    ProgressLogger progresslogger;
    progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
    progresslogger.startProgress(0, 1, "determining sequences to be removed");
    set<String> seq_to_remove;
    boost::random::mt19937 gen;
    while (seq_to_remove.size() < to_be_removed)
    {
      boost::random::uniform_int_distribution<> dist(0, maps[map_index].size()-1);
      int feature_idx = dist(gen);
      //OPENMS_LOG_INFO << "feature index chosen: " << feature_idx << std::endl;

      boost::random::uniform_int_distribution<> rand_pep_idx(0, maps[map_index][feature_idx].getPeptideIdentifications().size()-1);
      int pep_idx = rand_pep_idx(gen);
      //OPENMS_LOG_INFO << "peptide index chosen: " << pep_idx << std::endl;
      // get sequence from best hit
      seq_to_remove.insert(maps[map_index][feature_idx].getPeptideIdentifications()[pep_idx].getHits()[0].getSequence().toString());
    }
    progresslogger.endProgress();

    vector<FeatureMap> kept_maps(ins.size());
    vector<FeatureMap> removed_maps(ins.size());

    Size progress = 0;
    progresslogger.startProgress(0, ins.size(), "removing features from maps");
    for (Size i = 0; i < ins.size(); ++i)
    {
      progresslogger.setProgress(progress);
      // load map after map and clear after use because of copies
      if (i != map_index)
      {
        loadInputMap_(maps[i], ins[i], fxml_file);
      }
      //copy all properties
      kept_maps[i] = maps[i];
      //.. but delete feature information
      kept_maps[i].clear(false);
      removed_maps[i] = kept_maps[i];

      for (auto & feature : maps[i])
      {
        bool remove = false;
        for (vector<PeptideIdentification>::const_iterator pep_id_it = feature.getPeptideIdentifications().begin(); pep_id_it != feature.getPeptideIdentifications().end(); ++pep_id_it)
        {
          //loop over all peptideHits
          for (const auto & pep_hit_it : pep_id_it->getHits()) {
            if (seq_to_remove.find(pep_hit_it.getSequence().toString()) != seq_to_remove.end()) // PeptideHit contains sequence to be removed
            {
              remove = true;
            }
          }
        }
        // transfer feature to corresponding map
        if (remove)
        {
          removed_maps[i].push_back(feature);
        }
        else{
          kept_maps[i].push_back(feature);
        }
      }
      maps[i].clear(true);
      removed_maps[i].updateRanges();
      kept_maps[i].updateRanges();
      OPENMS_LOG_INFO << removed_maps[i].size() << " features from map " << (i+1) << " removed and " << kept_maps[i].size() << " features kept." << endl;
      ++progress;
    }
    progresslogger.endProgress();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    storeOutFiles_(kept_maps, outs, fxml_file);
    storeRemovedOutFiles_(removed_maps, removed_outs, fxml_file);

    return EXECUTION_OK;
  }

};

void TOPPSequenceRemover::loadInputMap_(FeatureMap& map, const String& in, FeatureXMLFile& fxml_file)
{
  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
  progresslogger.startProgress(0, 1, "loading input file "+in);
    fxml_file.load(in, map);
  progresslogger.endProgress();
}

void TOPPSequenceRemover::storeOutFiles_(const vector<FeatureMap>& out_maps, const StringList& out_files, FeatureXMLFile& fxml_file) {

  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
  progresslogger.startProgress(0, out_maps.size(), "writing output files");
  for (Size i = 0; i < out_files.size(); ++i)
  {
    progresslogger.setProgress(i);
    fxml_file.store(out_files[i], out_maps[i]);
  }
  progresslogger.endProgress();
}

void TOPPSequenceRemover::storeRemovedOutFiles_(const vector<FeatureMap>& removed_maps, const StringList& removed_out_files, FeatureXMLFile& fxml_file){
  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
  progresslogger.startProgress(0, removed_maps.size(), "writing files containing removed features");
  for (Size i = 0; i < removed_out_files.size(); ++i)
  {
    progresslogger.setProgress(i);
    fxml_file.store(removed_out_files[i], removed_maps[i]);
  }
  progresslogger.endProgress();
}


int main(int argc, const char** argv)
{
  TOPPSequenceRemover tool;
  return tool.main(argc, argv);
}

/// @endcond

