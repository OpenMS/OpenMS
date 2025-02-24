// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/FEATUREFINDER/SeedListGenerator.h>

#include <map>

// TODO REMOVE
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/SYSTEM/File.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_SeedListGenerator SeedListGenerator

@brief Application to generate seed lists for feature detection.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> &rarr; SeedListGenerator &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=3> @ref TOPP_FeatureFinderCentroided</td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping algorithm) </td>
        </tr>
    </table>
</CENTER>

Reference:\n
Weisser <em>et al.</em>: <a href="https://doi.org/10.1021/pr300992u">An automated pipeline for high-throughput label-free quantitative proteomics</a> (J. Proteome Res., 2013, PMID: 23391308).

In feature detection algorithms, an early step is generally to identify points of interest in the LC-MS map (so-called seeds) that may later be extended to features. If supported by the feature detection algorithm (currently only the "centroided" algorithm), user-supplied seed lists allow greater control over this process.

The SeedListGenerator can automatically create seed lists from a variety of sources. The lists are exported in featureXML format - suitable as input to FeatureFinder -, but can be converted to or from text formats using the @ref TOPP_TextExporter (with "-minimal" option to convert to CSV) and @ref TOPP_FileConverter (to convert from CSV) tools.


Seed lists can be generated from the file types below. The seeds are created at the indicated positions (RT/MZ):
<ul>
<li>mzML: locations of MS2 precursors
<li>idXML: locations of peptide identifications
<li>featureXML: locations of unassigned peptide identifications
<li>consensusXML: locations of consensus features that do not contain sub-features from the respective map
</ul>

If input is consensusXML, one output file per constituent map is required (same order as in the consensusXML). Otherwise, exactly one output file.

What are possible use cases for custom seed lists?
- In analyses that can take into account only features with peptide annotations, it may be useful to focus directly on certain locations in the LC-MS map - on all MS2 precursors (mzML input), or on precursors whose fragment spectra could be matched to a peptide sequence (idXML input).
- When additional information becomes available during an analysis, one might want to perform a second, targeted round of feature detection on the experimental data. For example, once a feature map is annotated with peptide identifications, it is possible to go back to the LC-MS map and look for features near unassigned peptides, potentially with a lower score threshold (featureXML input).
- Similarly, when features from different experiments are aligned and grouped, the consensus map may reveal where features were missed in the initial detection round in some experiments. The locations of these "holes" in the consensus map can be compiled into seed lists for the individual experiments (consensusXML input). (Note that the resulting seed lists use the retention time scale of the consensus map, which might be different from the original time scales of the experiments if e.g. one of the MapAligner tools was used to perform retention time correction as part of the alignment process. In this case, the RT transformations from the alignment must be applied to the LC-MS maps prior to the seed list-based feature detection runs.)

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SeedListGenerator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SeedListGenerator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{

  class TOPPSeedListGenerator :
    public TOPPBase
  {
public:
    TOPPSeedListGenerator() :
      TOPPBase("SeedListGenerator", "Generates seed lists for feature detection.")
    {
    }

protected:

    void registerOptionsAndFlags_() override
    {
      registerInputFile_("in", "<file>", "",
                         "Input file (see below for details)");
      setValidFormats_("in", ListUtils::create<String>("mzML,idXML,featureXML,consensusXML"));
      registerOutputPrefix_("out_prefix", "<prefix>", String(), "Output file prefix");
      setValidFormats_("out_prefix", ListUtils::create<String>("featureXML"));
      addEmptyLine_();
      registerFlag_("use_peptide_mass", "[idXML input only] Use the monoisotopic mass of the best peptide hit for the m/z position (default: use precursor m/z)");
    }

    ExitCodes main_(int, const char **) override
    {
      String in = getStringOption_("in");
      String out_prefix = getStringOption_("out_prefix");

      SeedListGenerator seed_gen;
      // results (actually just one result, except for consensusXML input):
      std::map<UInt64, SeedListGenerator::SeedList> seed_lists;

      Size num_maps = 0;
      FileTypes::Type in_type = FileHandler::getType(in);

      StringList out;
      out.push_back(out_prefix + "_0.featureXML"); // we manually set the name here

      if (in_type == FileTypes::CONSENSUSXML)
      {
        ConsensusMap consensus;
        FileHandler().loadConsensusFeatures(in, consensus, {FileTypes::CONSENSUSXML});
        num_maps = consensus.getColumnHeaders().size();
        ConsensusMap::ColumnHeaders ch = consensus.getColumnHeaders();
        size_t map_count = 0;
        // we have multiple out files
        out.clear();
        for([[maybe_unused]] const auto& header : ch)
        {           
          out.push_back(out_prefix + "_" + String(map_count) + ".featureXML"); // we manually set the name here
          ++map_count;
        }

        if (out.size() != num_maps)
        {
          writeLogError_("Error: expected " + String(num_maps) +
                    " output filenames");
          return ILLEGAL_PARAMETERS;
        }
        seed_gen.generateSeedLists(consensus, seed_lists);
      }
      else if (out.size() > 1)
      {
        writeLogError_("Error: expected only one output filename");
        return ILLEGAL_PARAMETERS;
      }
      else if (in_type == FileTypes::MZML)
      {
        PeakMap experiment;
        FileHandler().loadExperiment(in, experiment, {FileTypes::MZML});
        seed_gen.generateSeedList(experiment, seed_lists[0]);
      }
      else if (in_type == FileTypes::IDXML)
      {
        vector<ProteinIdentification> proteins;
        vector<PeptideIdentification> peptides;
        FileHandler().loadIdentifications(in, proteins, peptides, {FileTypes::IDXML});
        seed_gen.generateSeedList(peptides, seed_lists[0],
                                  getFlag_("use_peptide_mass"));
      }
      else if (in_type == FileTypes::FEATUREXML)
      {
        FeatureMap features;
        FileHandler().loadFeatures(in, features, {FileTypes::FEATUREXML});
        seed_gen.generateSeedList(
          features.getUnassignedPeptideIdentifications(), seed_lists[0]);
      }

      // output:
      num_maps = 0;
      for (std::map<UInt64, SeedListGenerator::SeedList>::iterator it =
             seed_lists.begin(); it != seed_lists.end(); ++it, ++num_maps)
      {
        FeatureMap features;
        seed_gen.convertSeedList(it->second, features);
        //annotate output with data processing info:
        addDataProcessing_(features, getProcessingInfo_(
                             DataProcessing::DATA_PROCESSING));
        OPENMS_LOG_INFO << "Writing " << features.size() << " seeds to " << out[num_maps] << endl;
        FileHandler().storeFeatures(out[num_maps], features, {FileTypes::FEATUREXML});
      }

      return EXECUTION_OK;
    }
  };
}

int main(int argc, const char ** argv)
{
  TOPPSeedListGenerator t;
  return t.main(argc, argv);
}

/// @endcond
