// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/AnnotatedMSRun.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

#include <type_traits>

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

  void removeDuplicates_(PeptideIdentificationList & peptides)
  {
    // there is no "PeptideIdentification::operator<", so we can't use a set
    // or sort + unique to filter out duplicates...
    // just use the naive O(n²) algorithm
    PeptideIdentificationList unique;
    for (PeptideIdentificationList::iterator in_it = peptides.begin();
         in_it != peptides.end(); ++in_it)
    {
      bool duplicate = false;
      for (PeptideIdentificationList::iterator out_it = unique.begin();
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

  /// Take the identifications of @p map (as peptide and protein identifications) and leave it without them
  template <typename MapType>
  static void split_(MapType& map, vector<ProteinIdentification>& proteins, PeptideIdentificationList& peptides)
  {
    IdentificationDataConverter::moveToIdentificationData(map);
    auto exported = IdentificationDataAdapter::toLegacy(map.getIdentificationData());
    proteins = std::move(exported.proteins);
    peptides = std::move(exported.peptides);
    const auto unlink = [](const auto& self, auto& feature) -> void {
      feature.getIDQueries().clear();
      feature.getIDMatches().clear();
      if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
      {
        for (auto& subordinate : feature.getSubordinates()) { self(self, subordinate); }
      }
    };
    for (auto& feature : map)
    {
      unlink(unlink, feature);
    }
    map.getIdentificationData().clear();
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (data annotated with identifications)");
    setValidFormats_("in", ListUtils::create<std::string>("featureXML,consensusXML"));
    registerOutputFile_("out", "<file>", "", "Output file (data without identifications). Either 'out' or 'id_out' are required. They can be used together.", false);
    setValidFormats_("out", ListUtils::create<std::string>("featureXML,consensusXML"));
    registerOutputFile_("id_out", "<file>", "", "Output file (identifications). Either 'out' or 'id_out' are required. They can be used together.", false);
    setValidFormats_("id_out", ListUtils::create<std::string>("idXML"));
  }

  ExitCodes main_(int, const char **) override
  {
    std::string in = getStringOption_("in"), out = getStringOption_("out"),
           id_out = getStringOption_("id_out");

    if (out.empty() && id_out.empty())
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__,
                                                 OPENMS_PRETTY_FUNCTION,
                                                 "out/id_out");
    }

    vector<ProteinIdentification> proteins;
    PeptideIdentificationList peptides;

    FileTypes::Type in_type = FileHandler::getType(in);

    if (in_type == FileTypes::FEATUREXML)
    {
      FeatureMap features;
      FileHandler().loadFeatures(in, features, {FileTypes::FEATUREXML});
      split_(features, proteins, peptides);
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
      split_(consensus, proteins, peptides);
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
      // resulting in duplicates
      removeDuplicates_(peptides);
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
