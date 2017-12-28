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
// $Maintainer: Timo Sachsenberg $
// $Authors: Alexandra Zerck, Chris Bielow$
// --------------------------------------------------------------------------

//#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/ANALYSIS/TARGETED/InclusionExclusionList.h>
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_InclusionExclusionListCreator InclusionExclusionListCreator

   @brief A tool for creating inclusion and/or exclusion lists for LC-MS/MS.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ InclusionExclusionListCreator \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> - </td>
        </tr>
        <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
        </tr>
    </table>
    </CENTER>

    Currently this tool can create tab-delimited inclusion or exclusion lists (m/z, RT start, RT stop).
    The input can either be peptide identifications from previous runs, a feature map or a FASTA-file with proteins.
    Inclusion and exclusion charges can be specified for FASTA and idXML input. If no charges are specified in the case of peptide id input, only
    the charge state of the peptide id is in/excluded, otherwise all given charge states are entered to the list.

    InclusionExclusionListCreator has different strategies for inclusion list creation: 'FeatureBased_LP', 'ProteinBased_LP' and 'ALL'.
    In the ALL mode all features are put onto the list. The FeatureBased_LP,  which was designed for MALDI data, maximizes the number of features in the inclusion list
    given the constraints that for each RT fraction a maximal number of precursors is not exceeded and each feature is scheduled at most a fixed number of times.
    In this mode, the sum of normalized feature intensities is maximized so that for each feature high intensity RTs are favoured over lower intensity ones.
    The ProteinBased_LP uses RT and detectability prediction methods to predict features that are most likely to be identified by MS/MS.
    Both LP methods are described in more detail in a recent publication: Zerck et al.:
    <a href="http://www.biomedcentral.com/1471-2105/14/56">Optimal precursor ion selection for LC-MALDI MS/MS</a> (BMC Bioinformatics 2013).

    The RT window size can be specified in the RT section of the INI file, either as relative window
    with [rt-rel_rt_window_size*rt,rt+rel_rt_window_size*rt] or absolute window.

    The default is RT in minutes, but seconds can also be used (see INI file).

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_InclusionExclusionListCreator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_InclusionExclusionListCreator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPInclusionExclusionListCreator :
  public TOPPBase
{
public:
  TOPPInclusionExclusionListCreator() :
    TOPPBase("InclusionExclusionListCreator", "Creates inclusion and/or exclusion lists.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("include", "<file>", "", "Inclusion list input file in FASTA or featureXML format.", false);
    setValidFormats_("include", ListUtils::create<String>("featureXML,fasta"));
    registerInputFile_("exclude", "<file>", "", "Exclusion list input file in featureXML, idXML or FASTA format.", false);
    setValidFormats_("exclude", ListUtils::create<String>("featureXML,idXML,fasta"));
    registerOutputFile_("out", "<file>", "", "Output file (tab delimited csv file).");
    setValidFormats_("out", ListUtils::create<String>("csv"));
    registerInputFile_("rt_model", "<file>", "", "RTModel file used for the rt prediction of peptides in FASTA files.", false);
    setValidFormats_("rt_model", ListUtils::create<String>("txt"));

    registerInputFile_("pt_model", "<file>", "", "PTModel file used for the pt prediction of peptides in FASTA files (only needed for inclusion_strategy PreotinBased_LP).", false);
    setValidFormats_("pt_model", ListUtils::create<String>("txt"));

    //in FASTA or featureXML
    registerIntList_("inclusion_charges", "<charge>", IntList(), "List containing the charge states to be considered for the inclusion list compounds, space separated.", false);
    setMinInt_("inclusion_charges", 1);
    registerStringOption_("inclusion_strategy", "<name>", "ALL", "strategy to be used for selection", false);
    setValidStrings_("inclusion_strategy", ListUtils::create<String>("FeatureBased_LP,ProteinBased_LP,ALL"));
    registerIntList_("exclusion_charges", "<charge>", IntList(), "List containing the charge states to be considered for the exclusion list compounds (for idXML and FASTA input), space separated.", false);
    setMinInt_("exclusion_charges", 1);
    registerInputFile_("raw_data", "<mzMLFile>", "", "File containing the raw data (only needed for FeatureBased_LP).", false);
    setValidFormats_("raw_data", ListUtils::create<String>("mzML"));

    //    setValidFormats_("out", ListUtils::create<String>("traML"));

    registerSubsection_("algorithm", "Inclusion/Exclusion algorithm section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    // there is only one subsection: 'algorithm' (s.a) .. and in it belongs the InclusionExclusionList param
    InclusionExclusionList fdc;
    OfflinePrecursorIonSelection ops;
    PSLPFormulation lp;
    Param tmp;
    tmp.insert("InclusionExclusionList:", fdc.getParameters());
    tmp.insert("PrecursorSelection:", ops.getParameters());
    tmp.remove("PrecursorSelection:mz_isolation_window");
    tmp.remove("PrecursorSelection:min_mz_peak_distance");
    tmp.insert("PrecursorSelection:", lp.getParameters().copy("feature_based"));
    return tmp;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input/output files
    String include(getStringOption_("include"));
    String exclude(getStringOption_("exclude"));
    String out(getStringOption_("out"));
    String strategy = getStringOption_("inclusion_strategy");
    std::cout << "strategy " << strategy << std::endl;
    String pt_model_file(getStringOption_("pt_model"));

    if (include == "" && exclude == "")
    {
      writeLog_("Error: No input file given.");
      return MISSING_PARAMETERS;
    }
    // currently we can handle only inclusion OR exclusion, will be possible with the traML output
    if (include != "" && exclude != "")
    {
      writeLog_("Error: Currently only inclusion OR exclusion, both will be possible with the traML output coming soon");
      return ILLEGAL_PARAMETERS;
    }

    IntList incl_charges(getIntList_("inclusion_charges"));
    IntList excl_charges(getIntList_("exclusion_charges"));
    String rt_model_file(getStringOption_("rt_model"));

    //-------------------------------------------------------------
    // loading input: inclusion list part
    //-------------------------------------------------------------

    FileHandler fh;
    TargetedExperiment exp;
    Param const& iel_param = getParam_().copy("algorithm:InclusionExclusionList:", true);
    writeDebug_("Parameters passed to InclusionExclusionList", iel_param, 3);

    InclusionExclusionList list;
    list.setParameters(iel_param);


    //    std::cout << "\n\n\n\n" << iel_param.getValue("RT:unit") << "\n\n";


    if (include != "")
    {
      FileTypes::Type in_type = fh.getType(include);

      if (in_type == FileTypes::FEATUREXML)
      {
        // load feature map
        FeatureMap map;
        FeatureXMLFile().load(include, map);

        if (strategy == "ALL")
        {
          if (!incl_charges.empty())
          {
            writeLog_("Warning: 'inclusion_charges' parameter is not honored for featureXML input with strategy ALL.");
            return ILLEGAL_PARAMETERS;
          }

          // convert to targeted experiment
          // for traML output
          //            list.loadTargets(map,incl_targets,exp);
          // for tab-delimited output
          try
          {
            list.writeTargets(map, out);
          }
          catch (Exception::UnableToCreateFile)
          {
            writeLog_("Error: Unable to create output file.");
            return CANNOT_WRITE_OUTPUT_FILE;
          }
        }
        else if (strategy == "FeatureBased_LP")
        {

          String raw_data_path = getStringOption_("raw_data");
          PeakMap exp, ms2;
          MzMLFile().load(raw_data_path, exp);
          FeatureMap out_map;
          StringList ms_runs;
          exp.getPrimaryMSRunPath(ms_runs);
          out_map.setPrimaryMSRunPath(ms_runs);

          IntList levels;
          levels.push_back(1);
          exp.getSpectra().erase(remove_if(exp.begin(), exp.end(),
                                           InMSLevelRange<MSSpectrum>(levels, true)), exp.end());
          exp.sortSpectra(true);
          OfflinePrecursorIonSelection opis;
          Param param = getParam_().copy("algorithm:PrecursorSelection:", true);
          param.removeAll("feature_based:");
          UInt spot_cap = param.getValue("ms2_spectra_per_rt_bin");
          opis.setParameters(param);

          // insert charges
          std::set<Int> charges_set;
          for (Size c = 0; c < incl_charges.size(); ++c)
          {
            charges_set.insert(incl_charges[c]);
          }

          // create ILP
          PSLPFormulation ilp_wrapper;
          Param param2 = getParam_().copy("algorithm:PrecursorSelection:", true);
          ilp_wrapper.setParameters(param2.copy("feature_based"));
          // get the mass ranges for each features for each scan it occurs in
          std::vector<std::vector<std::pair<Size, Size> > >  indices;
          opis.getMassRanges(map, exp, indices);

          std::vector<PSLPFormulation::IndexTriple> variable_indices;
          std::vector<int> solution_indices;
          ilp_wrapper.createAndSolveILPForKnownLCMSMapFeatureBased(map, exp, variable_indices, indices, charges_set, spot_cap, solution_indices);

          sort(variable_indices.begin(), variable_indices.end(), PSLPFormulation::IndexLess());
#ifdef DEBUG_OPS
          std::cout << "best_solution " << std::endl;
#endif
          std::vector<Int> rt_sizes(exp.size(), 0);
          // print best solution
          // create inclusion list
          for (Size i = 0; i < solution_indices.size(); ++i)
          {
            Size feature_index = variable_indices[solution_indices[i]].feature;
            Size scan = variable_indices[solution_indices[i]].scan;
            out_map.push_back(map[feature_index]);
            //            std::cout << map[feature_index].getMetaValue("msms_score") << std::endl;
            ++rt_sizes[scan];
          }
#ifdef DEBUG_OPS
          for (Size r = 0; r < rt_sizes.size(); ++r)
          {
            std::cout << r << "\t" << rt_sizes[r] << "\n";
          }
#endif
          try
          {
            if (out.hasSuffix("featureXML"))
            {
              FeatureXMLFile().store(out, out_map);
            }
            else list.writeTargets(out_map, out);
          }
          catch (Exception::UnableToCreateFile)
          {
            writeLog_("Error: Unable to create output file.");
            return CANNOT_WRITE_OUTPUT_FILE;
          }

        } //else if(strategy == "ILP") // ILP
        else
        {
          writeLog_("Warning: 'ProteinBased_LP' inclusion strategy is not valid for featureXML input.");
          return ILLEGAL_PARAMETERS;
        }
      }
      else // FASTA format
      {
        if (!File::exists(rt_model_file))
        {
          writeLog_("Error: RT model file required for FASTA input to predict RT elution time.");
          return MISSING_PARAMETERS;
        }
        if (incl_charges.empty())
        {
          writeLog_("Error: Protein sequences for inclusion given, but no charge states specified.");
          return MISSING_PARAMETERS;
        }
        if (strategy == "ProteinBased_LP")
        {
          OfflinePrecursorIonSelection opis;
          Param param = getParam_().copy("algorithm:PrecursorSelection:", true);
          param.removeAll("feature_based:");
          opis.setParameters(param);

          FeatureMap precursors;
          opis.createProteinSequenceBasedLPInclusionList(include, rt_model_file, pt_model_file, precursors);
          if (out.hasSuffix("featureXML"))
          {
            FeatureXMLFile().store(out, precursors);
          }
          else list.writeTargets(precursors, out);

        }
        else
        {
          std::vector<FASTAFile::FASTAEntry> entries;
          // load fasta-file
          FASTAFile().load(include, entries);

          // convert to targeted experiment
          // if traML output
          //list.loadTargets(entries,incl_targets,exp,missed_cleavages);
          // if tab-delimited output
          try
          {
            list.writeTargets(entries, out, incl_charges, rt_model_file);
          }
          catch (Exception::UnableToCreateFile)
          {
            writeLog_("Error: Unable to create output file.");
            return CANNOT_WRITE_OUTPUT_FILE;
          }
        }
      }

      //        exp.setIncludeTargets(incl_targets);
    }
    //-------------------------------------------------------------
    // loading input: exclusion list part
    //-------------------------------------------------------------
    if (exclude != "")
    {
      FileTypes::Type ex_type = fh.getType(exclude);
      //        std::vector<IncludeExcludeTarget> excl_targets;
      if (ex_type == FileTypes::FEATUREXML)
      {
        if (!excl_charges.empty())
        {
          writeLog_("Warning: 'exclusion_charges' parameter is not honored for featureXML input.");
          return ILLEGAL_PARAMETERS;
        }

        // load feature map
        FeatureMap map;
        FeatureXMLFile().load(exclude, map);

        // convert to targeted experiment if traML output is selected
        //            list.loadTargets(map,excl_targets,exp);
        // else write tab-delimited file directly
        try
        {
          list.writeTargets(map, out);
        }
        catch (Exception::UnableToCreateFile)
        {
          writeLog_("Error: Unable to create output file.");
          return CANNOT_WRITE_OUTPUT_FILE;
        }
      }
      else if (ex_type == FileTypes::IDXML)
      {
        std::vector<PeptideIdentification> pep_ids;
        std::vector<ProteinIdentification> prot_ids;
        IdXMLFile().load(exclude, prot_ids, pep_ids);
        try
        {
          list.writeTargets(pep_ids, out, excl_charges);
        }
        catch (Exception::UnableToCreateFile)
        {
          writeLog_("Error: Unable to create output file.");
          return CANNOT_WRITE_OUTPUT_FILE;
        }
        catch (Exception::InvalidSize)
        {
          writeLog_("Error: Peptide identification contains several hits. Use IDFilter to filter for significant peptide hits.");
          return ILLEGAL_PARAMETERS;
        }
        catch (Exception::MissingInformation)
        {
          writeLog_("Error: Peptide identification contains no RT information.");
          return ILLEGAL_PARAMETERS;
        }
      }
      else // FASTA format ...
      {
        if (!File::exists(rt_model_file))
        {
          writeLog_("Error: RT model file required for FASTA input to predict RT elution time.");
          return MISSING_PARAMETERS;
        }
        if (excl_charges.empty())
        {
          writeLog_("Error: Protein sequences for exclusion given, but no charge states specified.");
          return MISSING_PARAMETERS;
        }
        std::vector<FASTAFile::FASTAEntry> entries;
        // load fasta-file
        FASTAFile().load(exclude, entries);
        // convert to targeted experiment for traML output
        //            list.loadTargets(entries,excl_targets,exp,missed_cleavages);
        // else for tab-delimited output
        try
        {
          list.writeTargets(entries, out, excl_charges, rt_model_file);
        }
        catch (Exception::UnableToCreateFile)
        {
          writeLog_("Error: Unable to create output file.");
          return CANNOT_WRITE_OUTPUT_FILE;
        }
      }
//      exp.setExcludeTargets(excl_targets);
    }
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------


    //TraMLFile().store(out, exp);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPInclusionExclusionListCreator tool;
  return tool.main(argc, argv);
}

/// @endcond
