// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/OMSFileLoad.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FeatureFinderMetaboIdent FeatureFinderMetaboIdent

@brief Detects features in MS1 data corresponding to small molecule identifications.

<CENTER>
 <table>
   <tr>
     <td ALIGN="center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
     <td VALIGN="middle" ROWSPAN=2> &rarr; FeatureFinderMetaboIdent &rarr;</td>
     <td ALIGN="center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
   </tr>
   <tr>
     <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes (optional) </td>
     <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_TextExporter</td>
   </tr>
 </table>
</CENTER>

This tool detects quantitative features in MS1 data for a list of targets, typically small molecule/metabolite identifications.
It uses algorithms for targeted data analysis from the OpenSWATH pipeline.

@note This tool is still experimental!

@see @ref TOPP_FeatureFinderIdentification - targeted feature detection based on peptide identifications.

<B>Input format</B>

Spectra are expected in centroided or profile mode. Only MS1 level spectra are considered for feature detection.

The targets to quantify have to be specified in a tab-separated text file that is passed via the @p id parameter.
This file has to start with the following header line, defining its columns:
<pre>
<TT>CompoundName    SumFormula    Mass    Charge    RetentionTime    RetentionTimeRange    IsoDistribution</TT>
</pre>

Every subsequent line defines a target.
(Except lines starting with "#", which are considered as comments and skipped.)
The following requirements apply:
- @p CompoundName: unique name for the target compound
- @p SumFormula: chemical sum formula (see @ref OpenMS::EmpiricalFormula), optional
- @p Mass: neutral mass; if zero calculated from @p Formula
- @p Charge: charge state, or comma-separated list of multiple charges
- @p RetentionTime: retention time (RT), or comma-separated list of multiple RTs
- @p RetentionTimeRange: RT window around @p RetentionTime for chromatogram extraction, either one value or one per @p RT entry; if zero parameter @p extract:rt_window is used
- @p IsoDistribution: comma-separated list of relative abundances of isotopologues (see @ref OpenMS::IsotopeDistribution); if zero calculated from @p Formula

In the simplest case, only @p CompoundName, @p SumFormula, @p Charge and @p RetentionTime need to be given, all other values may be zero.
Every combination of compound (mass), RT and charge defines one target for feature detection.

<B>Output format</B>

The main output (parameter @p out) is a featureXML file containing the detected features, with annotations in meta data entries.
This file can be visualized in TOPPView - perhaps most usefully as a layer on top of the LC-MS data that gave rise to it.
Compound annotations of features (@p Name entries from the @p id input) can be shown by clicking the "Show feature annotation" button in the tool bar and selecting "Label meta data".
Positions of targets for which no feature was detected can be shown by clicking the "Show unassigned peptide identifications" button and selecting "Show label meta data".

To export the data from the featureXML file to a tabular text file (CSV), use @ref TOPP_TextExporter with the options @p no_ids and <TT>feature:add_metavalues 0</TT> (to include all meta data annotations).
In the result, the information from the @p CompoundName, @p SumFormula, @p Charge and @p RetentionTime columns from the input will be in the @p label, @p sum_formula, @p charge and @p expected_rt columns, respectively.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FeatureFinderMetaboIdent.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FeatureFinderMetaboIdent.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPFeatureFinderMetaboIdent :
  public TOPPBase
{
public:
  TOPPFeatureFinderMetaboIdent() :
    TOPPBase("FeatureFinderMetaboIdent", "Detects features in MS1 data based on metabolite identifications.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file: LC-MS raw data");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("id", "<file>", "", "Input file: Metabolite identifications");
    setValidFormats_("id", ListUtils::create<String>("tsv"));
    registerOutputFile_("out", "<file>", "", "Output file: Features");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerOutputFile_("lib_out", "<file>", "", "Output file: Assay library", false);
    setValidFormats_("lib_out", ListUtils::create<String>("traML"));
    registerOutputFile_("chrom_out", "<file>", "", "Output file: Chromatograms", false);
    setValidFormats_("chrom_out", ListUtils::create<String>("mzML"));
    registerOutputFile_("trafo_out", "<file>", "", "Output file: Retention times (expected vs. observed)", false);
    setValidFormats_("trafo_out", ListUtils::create<String>("trafoXML"));
    registerFlag_("force", "Force processing of files with no MS1 spectra", true);

    Param ffmetaboident_params;
    ffmetaboident_params.insert("", FeatureFinderAlgorithmMetaboIdent().getParameters());
    registerFullParam_(ffmetaboident_params); // register algorithm paramters as command line parameters
  }

  ProgressLogger prog_log_; ///< progress logger

  /// Read input file with information about targets
  vector<FeatureFinderAlgorithmMetaboIdent::FeatureFinderMetaboIdentCompound> readTargets_(const String& in_path)
  {
    vector<FeatureFinderAlgorithmMetaboIdent::FeatureFinderMetaboIdentCompound> metaboIdentTable;

    const string header =
      "CompoundName\tSumFormula\tMass\tCharge\tRetentionTime\tRetentionTimeRange\tIsoDistribution";
    ifstream source(in_path.c_str());
    if (!source.is_open())
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, in_path);
    }
    string line;
    getline(source, line);
    if (!String(line).hasPrefix(header))
    {
      String msg = "expected header line starting with: '" + header + "'";
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  line, msg);
    }
    Size line_count = 1;
    set<String> names;
    while (getline(source, line))
    {
      line_count++;
      if (line[0] == '#') continue; // skip comments
      vector<String> parts = ListUtils::create<String>(line, '\t'); // split
      if (parts.size() < 7)
      {
        OPENMS_LOG_ERROR
          << "Error: Expected 7 tab-separated fields, found only "
          << parts.size() << " in line " << line_count
          << " - skipping this line." << endl;
        continue;
      }
      String name = parts[0];
      if (name.empty())
      {
        OPENMS_LOG_ERROR << "Error: Empty name field in input line "
                         << line_count << " - skipping this line." << endl;
        continue;
      }
      if (!names.insert(name).second) // @TODO: is this check necessary?
      {
        OPENMS_LOG_ERROR << "Error: Duplicate name '" << name
                         << "' in input line " << line_count
                         << " - skipping this line." << endl;
        continue;
      }
      metaboIdentTable.push_back(FeatureFinderAlgorithmMetaboIdent::FeatureFinderMetaboIdentCompound(name,
                                 parts[1],
                                 parts[2].toDouble(),
                                 ListUtils::create<Int>(parts[3]),
                                 ListUtils::create<double>(parts[4]),
                                 ListUtils::create<double>(parts[5]),
                                 ListUtils::create<double>(parts[6])));
    }
    return metaboIdentTable;
  }

  /// Main function
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String id = getStringOption_("id");
    String out = getStringOption_("out");
    String lib_out = getStringOption_("lib_out");
    String chrom_out = getStringOption_("chrom_out");
    String trafo_out = getStringOption_("trafo_out");
    bool force = getFlag_("force");

    prog_log_.setLogType(log_type_);

    //-------------------------------------------------------------
    // load input
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Loading targets and creating assay library..." << endl;
    auto table = readTargets_(id);

    FeatureFinderAlgorithmMetaboIdent ff_mident;
    // copy (command line) tool parameters that match the algorithm parameters back to the algorithm
    auto tool_parameter = getParam_().copySubset(FeatureFinderAlgorithmMetaboIdent().getDefaults());
    tool_parameter.setValue("EMGScoring:init_mom", "true"); // overwrite defaults
    tool_parameter.setValue("EMGScoring:max_iteration", 100); // overwrite defaults
    tool_parameter.setValue("debug", debug_level_); // pass down debug level
    ff_mident.setParameters(tool_parameter);

    OPENMS_LOG_INFO << "Loading input LC-MS data..." << endl;
    FileHandler mzml;
    mzml.getOptions().addMSLevel(1);
    mzml.loadExperiment(in, ff_mident.getMSData(), {FileTypes::MZML});
    if (ff_mident.getMSData().empty() && !force)
    {
      OPENMS_LOG_ERROR << "Error: No MS1 scans in '"
                       << in << "' - aborting." << endl;
      return INCOMPATIBLE_INPUT_DATA;
    }
    FeatureMap features;
    ff_mident.run(table, features, in);

    // annotate "spectra_data" metavalue
    if (getFlag_("test"))
    {
      // if test mode set, add file without path so we can compare it
      features.setPrimaryMSRunPath({"file://" + File::basename(in)});
    }
    else
    {
      features.setPrimaryMSRunPath({in}, ff_mident.getMSData());
    }

    addDataProcessing_(features, getProcessingInfo_(DataProcessing::QUANTITATION));

    if (!chrom_out.empty())
    {
      PeakMap& chrom_data = ff_mident.getChromatograms();
      addDataProcessing_(chrom_data,
                         getProcessingInfo_(DataProcessing::FILTERING));
      FileHandler().storeExperiment(chrom_out, ff_mident.getChromatograms(), {FileTypes::MZML});
    }
    ff_mident.getChromatograms().clear(true);

    //-------------------------------------------------------------
    // write output
    //-------------------------------------------------------------

    OPENMS_LOG_INFO << "Writing final results..." << endl;
    FileHandler().storeFeatures(out, features, {FileTypes::FEATUREXML});

    // write transition library in TraML format
    if (!lib_out.empty())
    {
      FileHandler().storeTransitions(lib_out, ff_mident.getLibrary(), {FileTypes::TRAML});
    }

    // write expected vs. observed retention times
    if (!trafo_out.empty())
    {
      const TransformationDescription& trafo = ff_mident.getTransformations();
      FileHandler().storeTransformations(trafo_out, trafo, {FileTypes::TRANSFORMATIONXML});
    }

    //-------------------------------------------------------------
    // statistics
    //-------------------------------------------------------------

    Size n_missing = features.getUnassignedPeptideIdentifications().size();
    OPENMS_LOG_INFO << "\nSummary statistics:\n"
             << ff_mident.getLibrary().getCompounds().size() << " targets specified\n"
             << features.size() << " features found\n"
             << ff_mident.getNShared() << " features with multiple target annotations\n"
             << n_missing << " targets without features";
    const Size n_examples = 5;
    if (n_missing)
    {
      OPENMS_LOG_INFO << ":";
      for (Size i = 0;
           ((i < features.getUnassignedPeptideIdentifications().size()) &&
            (i < n_examples)); ++i)
      {
        const PeptideIdentification& id =
          features.getUnassignedPeptideIdentifications()[i];
        const TargetedExperiment::Compound& compound =
          ff_mident.getLibrary().getCompoundByRef(id.getMetaValue("PeptideRef"));
        OPENMS_LOG_INFO << "\n- " << ff_mident.prettyPrintCompound(compound);
      }
      if (n_missing > n_examples)
      {
        OPENMS_LOG_INFO << "\n- ... (" << n_missing - n_examples << " more)";
      }
    }
    OPENMS_LOG_INFO << "\n" << endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderMetaboIdent tool;
  return tool.main(argc, argv);
}

/// @endcond
