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

    @brief Removes user specified percentage of peptide identifications from a given featureXML file.

  Input and output format are 'featureXML'. The tool allows you to remove sequences
  from peptide identifications in featureXMLs. The removed sequences and their corresponding feature id will be
  written to the screen and can be parsed to a csv file.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSequenceRemover :
        public TOPPBase
{
public:
  TOPPSequenceRemover() :
          TOPPBase("SequenceRemover", "Removes user specified percentage of peptide identifications from each given featureXML file.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file is FeatureXML.");
    setValidFormats_("in", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out", "<file>", "", "Output file containing all features from input file, but lacks a user-specified percentage of peptide identifications.", false);
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerDoubleOption_("percent_to_remove", "<double>", 0.1, "Percentage of peptide identifications to be remove", false, false);
    setMinFloat_("percent_to_remove", 0.0);
    setMaxFloat_("percent_to_remove", 1.0);
  }

  ExitCodes main_(int, const char**) override
  {
    FeatureXMLFile fxml_file;
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    double percent_to_remove = getDoubleOption_("percent_to_remove");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    FeatureMap map;
    fxml_file.load(in, map);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    // determine amount of peptide ids to be removed
    boost::random::mt19937 gen;
    boost::random::uniform_int_distribution<> dist(0, map.size()-1);
    for (Size j = 0; j < (Size)ceil(map.size()*percent_to_remove); ++j)
    {
      int feature_idx = dist(gen);
      String feature_id = map[feature_idx].getUniqueId();
      Feature* f_ptr = &map[feature_idx];
      // store all sequences from PeptideIdentifications
      for (auto & p_it : f_ptr->getPeptideIdentifications())
      {
        cout << feature_id << '\t' << p_it.getHits()[0].getSequence().toString() << endl;
      }
      // remove PeptideIdentifications from Feature (all sequences)
      f_ptr->getPeptideIdentifications() = {};
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    fxml_file.store(out, map);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPSequenceRemover tool;
  return tool.main(argc, argv);
}

/// @endcond

