// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_IDSplitter IDSplitter

    @brief Splits protein/peptide identifications off of annotated data files.

    This performs the reverse operation as IDMapper.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_IDSplitter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_IDSplitter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDSplitter :
  public TOPPBase
{
public:

  TOPPIDSplitter() :
    TOPPBase("IDSplitter", "Splits protein/peptide identifications off of annotated data files", false)
  {
  }

protected:

  void removeDuplicates_(vector<PeptideIdentification> & peptides)
  {
    // there is no "PeptideIdentification::operator<", so we can't use a set
    // or sort + unique to filter out duplicates...
    // just use the naive O(n²) algorithm
    vector<PeptideIdentification> unique;
    for (vector<PeptideIdentification>::iterator in_it = peptides.begin();
         in_it != peptides.end(); ++in_it)
    {
      bool duplicate = false;
      for (vector<PeptideIdentification>::iterator out_it = unique.begin();
           out_it != unique.end(); ++out_it)
      {
        if (*in_it == *out_it)
        {
          duplicate = true;
          break;
        }
      }
      if (!duplicate) unique.push_back(*in_it);
    }
    peptides.swap(unique);
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (data annotated with identifications)");
    setValidFormats_("in", ListUtils::create<String>("mzML,featureXML,consensusXML"));
    registerOutputFile_("out", "<file>", "", "Output file (data without identifications). Either 'out' or 'id_out' are required. They can be used together.", false);
    setValidFormats_("out", ListUtils::create<String>("mzML,featureXML,consensusXML"));
    registerOutputFile_("id_out", "<file>", "", "Output file (identifications). Either 'out' or 'id_out' are required. They can be used together.", false);
    setValidFormats_("id_out", ListUtils::create<String>("idXML"));
  }

  ExitCodes main_(int, const char **) override
  {
    String in = getStringOption_("in"), out = getStringOption_("out"),
           id_out = getStringOption_("id_out");

    if (out.empty() && id_out.empty())
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__,
                                                 OPENMS_PRETTY_FUNCTION,
                                                 "out/id_out");
    }

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;

    FileTypes::Type in_type = FileHandler::getType(in);

    if (in_type == FileTypes::MZML)
    {
      PeakMap experiment;
      MzMLFile().load(in, experiment);
      // what about unassigned peptide IDs?
      for (PeakMap::Iterator exp_it = experiment.begin();
           exp_it != experiment.end(); ++exp_it)
      {
        peptides.insert(peptides.end(),
                        exp_it->getPeptideIdentifications().begin(),
                        exp_it->getPeptideIdentifications().end());
        exp_it->getPeptideIdentifications().clear();
      }
      experiment.getProteinIdentifications().swap(proteins);
      if (!out.empty())
      {
        addDataProcessing_(experiment,
                           getProcessingInfo_(DataProcessing::FILTERING));
        MzMLFile().store(out, experiment);
      }
    }
    else if (in_type == FileTypes::FEATUREXML)
    {
      FeatureMap features;
      FeatureXMLFile().load(in, features);
      features.getUnassignedPeptideIdentifications().swap(peptides);
      for (FeatureMap::Iterator feat_it = features.begin();
           feat_it != features.end(); ++feat_it)
      {
        peptides.insert(peptides.end(),
                        feat_it->getPeptideIdentifications().begin(),
                        feat_it->getPeptideIdentifications().end());
        feat_it->getPeptideIdentifications().clear();
      }
      features.getProteinIdentifications().swap(proteins);
      if (!out.empty())
      {
        addDataProcessing_(features,
                           getProcessingInfo_(DataProcessing::FILTERING));
        FeatureXMLFile().store(out, features);
      }
    }
    else         // consensusXML
    {
      ConsensusMap consensus;
      ConsensusXMLFile().load(in, consensus);
      consensus.getUnassignedPeptideIdentifications().swap(peptides);
      for (ConsensusMap::Iterator cons_it = consensus.begin();
           cons_it != consensus.end(); ++cons_it)
      {
        peptides.insert(peptides.end(),
                        cons_it->getPeptideIdentifications().begin(),
                        cons_it->getPeptideIdentifications().end());
        cons_it->getPeptideIdentifications().clear();
      }
      consensus.getProteinIdentifications().swap(proteins);
      if (!out.empty())
      {
        addDataProcessing_(consensus,
                           getProcessingInfo_(DataProcessing::FILTERING));
        ConsensusXMLFile().store(out, consensus);
      }
    }

    if (!id_out.empty())
    {
      // IDMapper can match a peptide ID to several overlapping features,
      // resulting in duplicates; this shouldn't be the case for peak data
      if (in_type != FileTypes::MZML) removeDuplicates_(peptides);
      IdXMLFile().store(id_out, proteins, peptides);
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPIDSplitter tool;
  return tool.main(argc, argv);
}

/// @endcond
