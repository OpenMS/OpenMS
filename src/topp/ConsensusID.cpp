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
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPMatrix.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPIons.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmBest.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmWorst.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmAverage.h>
#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_ConsensusID ConsensusID

    @brief Computes a consensus from results of multiple peptide identification engines.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ ConsensusID \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDPosteriorErrorProbability </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=3> @ref TOPP_PeptideIndexer </td>
        </tr>
        <tr>
          <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_IDFilter </td>
        </tr>
        <tr>
          <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_IDMapper </td>
        </tr>
    </table>
    </CENTER>

    <B>Reference:</B>

    Nahnsen <em>et al.</em>: <a href="http://dx.doi.org/10.1021/pr2002879">Probabilistic consensus scoring improves tandem mass spectrometry peptide identification</a> (J. Proteome Res., 2011, PMID: 21644507).

    <B>Algorithms:</B>

    ConsensusID offers several algorithms that can aggregate results from multiple peptide identification engines ("search engines") into consensus identifications - typically one per MS2 spectrum. This works especially well for search engines that provide more than one peptide hit per spectrum, i.e. that report not just the best hit, but also a list of runner-up candidates with corresponding scores.

    The available algorithms are (see also @ref OpenMS::ConsensusIDAlgorithm and its subclasses):
    @li @p PEPMatrix: Scoring based on posterior error probabilities (PEPs) and peptide sequence similarities. This algorithm uses a substitution matrix to score the similarity of sequences not listed by all search engines. It requires PEPs as the scores for all peptide hits.
    @li @p PEPIons: Scoring based on posterior error probabilities (PEPs) and fragment ion similarities ("shared peak count"). This algorithm, too, requires PEPs as scores.
    @li @p best: For each peptide ID, this uses the best score of any search engine as the consensus score. All peptide IDs must have the same score type.
    @li @p worst: For each peptide ID, this uses the worst score of any search engine as the consensus score. All peptide IDs must have the same score type.
    @li @p average: For each peptide ID, this uses the average score of all search engines as the consensus score. Again, all peptide IDs must have the same score type.
    @li @p ranks: Calculates a consensus score based on the ranks of peptide IDs in the results of different search engines. The final score is in the range (0, 1], with 1 being the best score. The input peptide IDs do not need to have the same score type.

    PEPs for search results can be calculated using the @ref TOPP_IDPosteriorErrorProbability tool, which supports a variety of search engines.

    @note Important: All protein-level identification results will be lost by applying ConsensusID. (It is unclear how potentially conflicting protein-level results from different search engines should be combined.) If necessary, run the @ref TOPP_PeptideIndexer tool to add protein references for peptides again.

    @note Peptides with different post-translational modifications (PTMs), or with different site localizations of the same PTMs, are treated as different peptides by all algorithms. However, a qualification applies for the @p PEPMatrix algorithm: The similarity scoring method used there can only take unmodified peptide sequences into account, so PTMs are ignored during that step. However, the PTMs are not removed from the peptides, and there will be separate results for differently-modified peptides.

    <B>File types:</B>

    Different input files types are supported:
    @li idXML: A file containing multiple identification runs, typically from different search engines. Use @ref TOPP_IDMerger to merge individual idXML files from different search runs into one. During the ConsensusID analysis, the identification results will be grouped according to their originating MS2 spectra, based on retention time and precursor m/z information (see parameters @p rt_delta and @p mz_delta). One consensus identification will be generated for each group.
    @li featureXML or consensusXML: Given (consensus) features annotated with peptide identifications from multiple search runs, one consensus identification is created for every annotated feature. Peptide identifications not assigned to features are not considered and will be removed. See @ref TOPP_IDMapper for the task of mapping peptide identifications to feature maps or consensus maps.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>Filtering:</B>

    Generally, search results can be filtered according to various criteria using @ref TOPP_IDFilter before (or after) applying this tool. ConsensusID itself offers only a limited number of filtering options that are especially useful in its context (see the @p filter parameter section):
    @li @p considered_hits: Limits the number of alternative peptide hits considered per spectrum/feature for each identification run. This helps to reduce runtime, especially for the @p PEPMatrix and @p PEPIons algorithms, which involve costly "all vs. all" comparisons of peptide hits.
    @li @p min_support: This allows filtering of peptide hits based on agreement between search engines. Every peptide sequence in the analysis has been identified by at least one search run. This parameter defines which fraction (between 0 and 1) of the remaining search runs must "support" a peptide identification that should be kept. The meaning of "support" differs slightly between algorithms: For @p best, @p worst, @p average and @p rank, each search run supports peptides that it has also identified among its top @p considered_hits candidates. So @p min_support simply gives the fraction of additional search engines that must have identified a peptide. (For example, if there are three search runs, and only peptides identified by at least two of them should be kept, set @p min_support to 0.5.) For the similarity-based algorithms @p PEPMatrix and @p PEPIons, the "support" for a peptide is the average similarity of the most-similar peptide from each (other) search run. (In the context of the JPR publication, this is the average of the similarity scores used in the consensus score calculation for a peptide.)
    @li @p count_empty: Typically not all search engines will provide results for all searched MS2 spectra. This parameter determines whether search runs that provided no results should be counted in the "support" calculation; by default, they are ignored.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_ConsensusID.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_ConsensusID.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPConsensusID :
  public TOPPBase
{
public:
  TOPPConsensusID() :
    TOPPBase("ConsensusID", "Computes a consensus of peptide identifications of several identification engines.")
  {
  }

protected:

  String algorithm_; // algorithm for consensus calculation (input parameter)

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("idXML,featureXML,consensusXML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("idXML,featureXML,consensusXML"));

    addEmptyLine_();
    registerDoubleOption_("rt_delta", "<value>", 0.1, "[idXML input only] Maximum allowed retention time deviation between identifications belonging to the same spectrum.", false);
    setMinFloat_("rt_delta", 0.0);
    registerDoubleOption_("mz_delta", "<value>", 0.1, "[idXML input only] Maximum allowed precursor m/z deviation between identifications belonging to the same spectrum.", false);
    setMinFloat_("mz_delta", 0.0);

    // General algorithm parameters are defined in the abstract base class
    // "ConsensusIDAlgorithm", but we can't get them from there because we can't
    // instantiate the class. So we get those parameters from a subclass that
    // doesn't add any other parameters:
    registerTOPPSubsection_("filter", "Options for filtering peptide hits");
    registerFullParam_(ConsensusIDAlgorithmBest().getDefaults());

    registerStringOption_("algorithm", "<choice>", "PEPMatrix",
                          "Algorithm used for consensus scoring.\n"
                          "* PEPMatrix: Scoring based on posterior error probabilities (PEPs) and peptide sequence similarities (scored by a substitution matrix). Requires PEPs as scores.\n"
                          "* PEPIons: Scoring based on posterior error probabilities (PEPs) and fragment ion similarities ('shared peak count'). Requires PEPs as scores.\n"
                          "* best: For each peptide ID, use the best score of any search engine as the consensus score. Requires the same score type in all ID runs.\n"
                          "* worst: For each peptide ID, use the worst score of any search engine as the consensus score. Requires the same score type in all ID runs.\n"
                          "* average:  For each peptide ID, use the average score of all search engines as the consensus. Requires the same score type in all ID runs.\n"
                          "* ranks: Calculates a consensus score based on the ranks of peptide IDs in the results of different search engines. The final score is in the range (0, 1], with 1 being the best score. No requirements about score types.", false);
    setValidStrings_("algorithm", ListUtils::create<String>("PEPMatrix,PEPIons,best,worst,average,ranks"));

    // subsections appear in alphabetical (?) order, independent of the order
    // in which they were registered:
    registerSubsection_("PEPIons", "PEPIons algorithm parameters");
    registerSubsection_("PEPMatrix", "PEPMatrix algorithm parameters");
  }


  Param getSubsectionDefaults_(const String& section) const override
  {
    Param algo_params;
    if (section == "PEPMatrix")
    {
      algo_params = ConsensusIDAlgorithmPEPMatrix().getDefaults();
    }
    else // section == "PEPIons"
    {
      algo_params = ConsensusIDAlgorithmPEPIons().getDefaults();
    }
    // remove parameters defined in the base class (to avoid duplicates):
    algo_params.remove("filter:");
    return algo_params;
  }


  void setProteinIdentifications_(vector<ProteinIdentification>& prot_ids)
  {
    // modification params are necessary for further analysis tools (e.g. LuciPHOr2)
    set<String> fixed_mods_set;
    set<String> var_mods_set;
    StringList merged_spectra_data;
    for (vector<ProteinIdentification>::iterator it_prot_ids = prot_ids.begin(); it_prot_ids != prot_ids.end(); ++it_prot_ids)
    {
      ProteinIdentification::SearchParameters search_params(it_prot_ids->getSearchParameters());
      std::copy(search_params.fixed_modifications.begin(), search_params.fixed_modifications.end(), std::inserter(fixed_mods_set, fixed_mods_set.end()));
      std::copy(search_params.variable_modifications.begin(), search_params.variable_modifications.end(), std::inserter(var_mods_set, var_mods_set.end()));
      StringList spectra_data;
      it_prot_ids->getPrimaryMSRunPath(spectra_data);
      std::copy(spectra_data.begin(), spectra_data.end(), std::inserter(merged_spectra_data, merged_spectra_data.end()));
    }
    ProteinIdentification::SearchParameters search_params;
    std::vector<String> fixed_mods(fixed_mods_set.begin(), fixed_mods_set.end());
    std::vector<String> var_mods(var_mods_set.begin(), var_mods_set.end());
    search_params.fixed_modifications    = fixed_mods;
    search_params.variable_modifications = var_mods;

    prot_ids.clear();
    prot_ids.resize(1);
    prot_ids[0].setDateTime(DateTime::now());
    prot_ids[0].setSearchEngine("OpenMS/ConsensusID_" + algorithm_);
    prot_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    prot_ids[0].setSearchParameters(search_params);

    // make file name entries unique
    std::sort(merged_spectra_data.begin(), merged_spectra_data.end());
    StringList::iterator last = std::unique(merged_spectra_data.begin(), merged_spectra_data.end());
    merged_spectra_data.erase(last, merged_spectra_data.end());
    prot_ids[0].setPrimaryMSRunPath(merged_spectra_data);
  }


  template <typename MapType>
  void processFeatureOrConsensusMap_(MapType& input_map,
                                     ConsensusIDAlgorithm* consensus)
  {
    // Problem with feature data: IDs from multiple spectra may be attached to
    // a (consensus) feature, so we may have multiple IDs from the same search
    // engine. This means that we can't just use the number of search runs as
    // our "baseline" for the number of identifications (parameter
    // "number_of_runs")! To work around this, we multiply the number of
    // different ID runs with the max. number of times we see the same ID run
    // in the annotations of a feature.

    map<String, Size> id_mapping; // mapping: run ID -> index
    Size number_of_runs = input_map.getProteinIdentifications().size();
    for (Size i = 0; i < number_of_runs; ++i)
    {
      id_mapping[input_map.getProteinIdentifications()[i].getIdentifier()] = i;
    }

    // compute consensus:
    for (typename MapType::Iterator map_it = input_map.begin();
         map_it != input_map.end(); ++map_it)
    {
      vector<PeptideIdentification>& ids = map_it->getPeptideIdentifications();
      vector<Size> times_seen(number_of_runs);
      for (vector<PeptideIdentification>::iterator pep_it = ids.begin(); 
           pep_it != ids.end(); ++pep_it)
      {
        ++times_seen[id_mapping[pep_it->getIdentifier()]];
      }
      Size n_repeats = *max_element(times_seen.begin(), times_seen.end());

      consensus->apply(ids, number_of_runs * n_repeats);
    }

    // create new identification run:
    setProteinIdentifications_(input_map.getProteinIdentifications());
    // remove outdated information (protein references will be broken):
    input_map.getUnassignedPeptideIdentifications().clear();
  }


  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in");
    FileTypes::Type in_type = FileHandler::getType(in);
    String out = getStringOption_("out");
    double rt_delta = getDoubleOption_("rt_delta");
    double mz_delta = getDoubleOption_("mz_delta");

    //----------------------------------------------------------------
    // set up ConsensusID
    //----------------------------------------------------------------
    ConsensusIDAlgorithm* consensus;
    // general algorithm parameters:
    Param algo_params = ConsensusIDAlgorithmBest().getDefaults();
    algorithm_ = getStringOption_("algorithm");
    if (algorithm_ == "PEPMatrix")
    {
      consensus = new ConsensusIDAlgorithmPEPMatrix();
      // add algorithm-specific parameters:
      algo_params.merge(getParam_().copy("PEPMatrix:", true));
    }
    else if (algorithm_ == "PEPIons")
    {
      consensus = new ConsensusIDAlgorithmPEPIons();
      // add algorithm-specific parameters:
      algo_params.merge(getParam_().copy("PEPIons:", true));
    }
    else if (algorithm_ == "best")
    {
      consensus = new ConsensusIDAlgorithmBest();
    }
    else if (algorithm_ == "worst")
    {
      consensus = new ConsensusIDAlgorithmWorst();
    }
    else if (algorithm_ == "average")
    {
      consensus = new ConsensusIDAlgorithmAverage();
    }
    else // algorithm_ == "ranks"
    {
      consensus = new ConsensusIDAlgorithmRanks();
    }
    algo_params.update(getParam_(), false, Log_debug); // update general params.
    consensus->setParameters(algo_params);

    //----------------------------------------------------------------
    // idXML
    //----------------------------------------------------------------
    if (in_type == FileTypes::IDXML)
    {
      vector<ProteinIdentification> prot_ids;
      vector<PeptideIdentification> pep_ids;
      String document_id;
      IdXMLFile().load(in, prot_ids, pep_ids, document_id);

      // merge peptide IDs by precursor position - this is equivalent to a
      // feature linking problem (peptide IDs from different ID runs <-> 
      // features from different maps), so we bring the data into a format
      // suitable for a feature grouping algorithm:
      vector<FeatureMap> maps(prot_ids.size());
      map<String, Size> id_mapping; // mapping: run ID -> index (of feature map)
      for (Size i = 0; i < prot_ids.size(); ++i)
      {
        id_mapping[prot_ids[i].getIdentifier()] = i;
      }

      for (vector<PeptideIdentification>::iterator pep_it = pep_ids.begin();
           pep_it != pep_ids.end(); ++pep_it)
      {
        String run_id = pep_it->getIdentifier();
        if (!pep_it->hasRT() || !pep_it->hasMZ())
        {
          LOG_FATAL_ERROR << "Peptide ID without RT and/or m/z information found in identification run '" + run_id + "'.\nMake sure that this information is included for all IDs when generating/converting search results. Aborting!" << endl;
          return INCOMPATIBLE_INPUT_DATA;
        }

        Feature feature;
        feature.setRT(pep_it->getRT());
        feature.setMZ(pep_it->getMZ());
        feature.getPeptideIdentifications().push_back(*pep_it);
        maps[id_mapping[run_id]].push_back(feature);
      }
      // precondition for "FeatureGroupingAlgorithmQT::group":
      for (vector<FeatureMap>::iterator map_it = maps.begin();
           map_it != maps.end(); ++map_it)
      {
        map_it->updateRanges();
      }

      FeatureGroupingAlgorithmQT linker;
      Param linker_params = linker.getDefaults();
      linker_params.setValue("use_identifications", "false");
      linker_params.setValue("ignore_charge", "true");
      linker_params.setValue("distance_RT:max_difference", rt_delta);
      linker_params.setValue("distance_MZ:max_difference", mz_delta);
      linker_params.setValue("distance_MZ:unit", "Da");
      linker.setParameters(linker_params);

      ConsensusMap grouping;
      linker.group(maps, grouping);

      // compute consensus
      pep_ids.clear();
      for (ConsensusMap::Iterator it = grouping.begin(); it != grouping.end();
           ++it)
      {
        consensus->apply(it->getPeptideIdentifications(), prot_ids.size());
        if (!it->getPeptideIdentifications().empty())
        {
          PeptideIdentification& pep_id = it->getPeptideIdentifications()[0];
          // hits may be empty due to filtering (parameter "min_support");
          // in that case skip to avoid a warning from "IDXMLFile::store":
          if (!pep_id.getHits().empty())
          {
            pep_id.setRT(it->getRT());
            pep_id.setMZ(it->getMZ());
            pep_ids.push_back(pep_id);
          }
        }
      }

      // create new identification run
      setProteinIdentifications_(prot_ids);

      // store consensus
      IdXMLFile().store(out, prot_ids, pep_ids);
    }


    //----------------------------------------------------------------
    // featureXML
    //----------------------------------------------------------------
    if (in_type == FileTypes::FEATUREXML)
    {
      FeatureMap map;
      FeatureXMLFile().load(in, map);

      processFeatureOrConsensusMap_(map, consensus);

      FeatureXMLFile().store(out, map);
    }

    //----------------------------------------------------------------
    // consensusXML
    //----------------------------------------------------------------
    if (in_type == FileTypes::CONSENSUSXML)
    {
      ConsensusMap map;
      ConsensusXMLFile().load(in, map);

      processFeatureOrConsensusMap_(map, consensus);

      ConsensusXMLFile().store(out, map);
    }

    delete consensus;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPConsensusID tool;
  return tool.main(argc, argv);
}

/// @endcond
