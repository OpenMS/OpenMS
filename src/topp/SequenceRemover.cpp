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
          TOPPMapAlignerBase("SequenceRemover", "Removes 'n' sequences from peptides randomly from featureXML files.")
  {

  }

  // function declarations
  void loadInputMaps_(vector<FeatureMap>& maps, StringList& ins, FeatureXMLFile& fxml_file);
  void storeFeatureXMLs_(const vector<FeatureMap>& feature_maps, const StringList& out_files, FeatureXMLFile& fxml_file);

protected:
  void registerOptionsAndFlags_() override
  {
    TOPPMapAlignerBase::registerOptionsAndFlags_("featureXML", REF_NONE);
    /*
    registerInputFile_("ins", "<files>", "", "input files");
    setValidFormats_("ins", {"featureXML"});
    registerOutputFile_("outs", "<files>", "", "output files");
    setValidFormats_("outs", {"featureXML"});
     */
    registerIntOption_("number_of_sequences", "<int>", 1, "Number of randomly chosen peptide sequences to remove", false);
    setMinInt_("number_of_sequences", 1);
  }

  ExitCodes main_(int, const char**) override
  {
    FeatureXMLFile fxml_file;
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList ins = getStringList_("in");
    StringList outs = getStringList_("out");
    Size number_of_sequences = getIntOption_("number_of_sequences");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    vector<FeatureMap> maps(ins.size());
    loadInputMaps_(maps, ins, fxml_file);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    vector<String> to_remove(number_of_sequences);
    for (Size i = 0; i < number_of_sequences; ++i)
    {
      Size rand = std::rand();
      std::cout << "rand int: " << rand << std::endl;
      Size rand_map_idx = rand % maps.size();
      Size rand_feature_idx = rand % maps[rand_map_idx].size();
      Size rand_pep_idx = rand % maps[rand_map_idx][rand_feature_idx].getPeptideIdentifications().size();
      // remove sequence from best hit
      to_remove[i] = maps[rand_map_idx][rand_feature_idx].getPeptideIdentifications()[rand_pep_idx].getHits()[0].getSequence().toString();
    }
    AASequence empty_seq;
    empty_seq.fromString("");

    for ( auto & map : maps )
    {
      for (auto & feature : map)
      {
        for (auto pep_it = feature.getPeptideIdentifications().begin(); pep_it != feature.getPeptideIdentifications().end(); ++pep_it)
        {
          if (! pep_it->getHits().empty())
          {
            for (auto hit_it = pep_it->getHits().begin(); hit_it != pep_it->getHits().end(); ++hit_it)
            {
              for (auto & seq : to_remove)
              {
                if (hit_it->getSequence().toString() == seq)
                {
                  hit_it->setSequence(empty_seq);
                  OPENMS_LOG_INFO << "removed sequence: " << seq << std::endl;
                }
              }
            }
          }
        }
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    storeFeatureXMLs_(maps, outs, fxml_file);

    return EXECUTION_OK;
  }

};

void TOPPSequenceRemover::loadInputMaps_(vector<FeatureMap>& maps, StringList& ins, FeatureXMLFile& fxml_file)
{
  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
  progresslogger.startProgress(0, ins.size(), "loading input files");
  for (Size i = 0; i < ins.size(); ++i)
  {
    progresslogger.setProgress(i);
    fxml_file.load(ins[i], maps[i]);
  }
  progresslogger.endProgress();
}

void TOPPSequenceRemover::storeFeatureXMLs_(const vector<FeatureMap>& feature_maps, const StringList& out_files, FeatureXMLFile& fxml_file) {

  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
  progresslogger.startProgress(0, feature_maps.size(), "writing output files");
  for (Size i = 0; i < out_files.size(); ++i)
  {
    progresslogger.setProgress(i);
    fxml_file.store(out_files[i], feature_maps[i]);
  }
  progresslogger.endProgress();
}


int main(int argc, const char** argv)
{
  TOPPSequenceRemover tool;
  return tool.main(argc, argv);
}

/// @endcond

