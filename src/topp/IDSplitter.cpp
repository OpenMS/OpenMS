// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IDSplitter IDSplitter

@brief Splits protein/peptide identifications off of annotated data files.

This performs the reverse operation as IDMapper.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IDSplitter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IDSplitter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDSplitter :
  public TOPPBase
{
public:

  TOPPIDSplitter() :
    TOPPBase("IDSplitter", "Splits protein/peptide identifications off of annotated data files")
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
      FileHandler().loadExperiment(in, experiment, {FileTypes::MZML});
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
        FileHandler().storeExperiment(out, experiment, {FileTypes::MZML});
      }
    }
    else if (in_type == FileTypes::FEATUREXML)
    {
      FeatureMap features;
      FileHandler().loadFeatures(in, features, {FileTypes::FEATUREXML});
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
        FileHandler().storeFeatures(out, features, {FileTypes::FEATUREXML});
      }
    }
    else         // consensusXML
    {
      ConsensusMap consensus;
      FileHandler().loadConsensusFeatures(in, consensus, {FileTypes::CONSENSUSXML});
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
        FileHandler().storeConsensusFeatures(out, consensus, {FileTypes::CONSENSUSXML});
      }
    }

    if (!id_out.empty())
    {
      // IDMapper can match a peptide ID to several overlapping features,
      // resulting in duplicates; this shouldn't be the case for peak data
      if (in_type != FileTypes::MZML) removeDuplicates_(peptides);
      FileHandler().storeIdentifications(id_out, proteins, peptides, {FileTypes::IDXML});
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
