// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
#include <cstdlib>  // for "rand"
#include <ctime>    // for "time" (seeding of random number generator)
#include <iterator> // for "inserter", "back_inserter"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FeatureFinderIdentification FeatureFinderIdentification

@brief Detects features in MS1 data based on peptide identifications.

<CENTER>
 <table>
   <tr>
     <th ALIGN = "center"> pot. predecessor tools </td>
     <td VALIGN="middle" ROWSPAN=3> &rarr; FeatureFinderIdentification &rarr;</td>
     <th ALIGN = "center"> pot. successor tools </td>
   </tr>
   <tr>
     <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes (optional) </td>
     <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_ProteinQuantifier</td>
   </tr>
   <tr>
     <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
   </tr>
 </table>
</CENTER>

@b Reference: @n
Weisser & Choudhary: <a href="https://doi.org/10.1021/acs.jproteome.7b00248">Targeted Feature Detection for Data-Dependent Shotgun Proteomics</a> (J. Proteome Res., 2017, PMID: 28673088).

This tool detects quantitative features in MS1 data based on information from peptide identifications (derived from MS2 spectra).
It uses algorithms for targeted data analysis from the OpenSWATH pipeline.

The aim is to detect features that enable the quantification of (ideally) all peptides in the identification input.
This is based on the following principle: When a high-confidence identification (ID) of a peptide was made based on an MS2 spectrum from a certain (precursor) position in the LC-MS map, this
indicates that the particular peptide is present at that position, so a feature for it should be detectable there.

@note It is important that only high-confidence (i.e. reliable) peptide identifications are used as input!

Targeted data analysis on the MS1 level uses OpenSWATH algorithms and follows roughly the steps outlined below.

<B>Use of inferred ("external") IDs</B>

The situation becomes more complicated when several LC-MS/MS runs from related samples of a label-free experiment are considered.
In order to quantify a larger fraction of the peptides/proteins in the samples, it is desirable to infer peptide identifications across runs.
Ideally, all peptides identified in any of the runs should be quantified in each and every run.
However, for feature detection of inferred ("external") IDs, the following problems arise:
First, retention times may be shifted between the run being quantified and the run that gave rise to the ID.
Such shifts can be corrected (see @ref TOPP_MapAlignerIdentification), but only to an extent.
Thus, the RT location of the inferred ID may not necessarily lie within the RT range of the correct feature.
Second, since the peptide in question was not directly identified in the run being quantified, it may not actually be present in detectable amounts in that sample, e.g. due to differential
regulation of the corresponding protein. There is thus a risk of introducing false-positive features.

FeatureFinderIdentification deals with these challenges by explicitly distinguishing between internal IDs (derived from the LC-MS/MS run being quantified) and external IDs (inferred from related
runs). Features derived from internal IDs give rise to a training dataset for an SVM classifier. The SVM is then used to predict which feature candidates derived from external IDs are most likely
to be correct. See steps 4 and 5 below for more details.

<B>1. Assay generation</B>

Feature detection is based on assays for identified peptides, each of which incorporates the retention time (RT), mass-to-charge ratio (m/z), and isotopic distribution (derived from the sequence)
of a peptide. Peptides with different modifications are considered different peptides. One assay will be generated for every combination of (modified) peptide sequence, charge state, and RT region
that has been identified. The RT regions arise by pooling all identifications of the same peptide, considering a window of size @p extract:rt_window around every RT location that gave rise to an
ID, and then merging overlapping windows.

<B>2. Ion chromatogram extraction</B>

Ion chromatograms (XICs) are extracted from the LC-MS data (parameter @p in).
One XIC per isotope in an assay is generated, with the corresponding m/z value and RT range (variable, depending on the RT region of the assay).

@see @ref TOPP_OpenSwathChromatogramExtractor

<B>3. Feature detection</B>

Next feature candidates - typically several per assay - are detected in the XICs and scored.
A variety of scores for different quality aspects are calculated by OpenSWATH.

@see @ref TOPP_OpenSwathAnalyzer

<B>4. Feature classification</B>

Feature candidates derived from assays with "internal" IDs are classed as "negative" (candidates without matching internal IDs), "positive" (the single best candidate per assay with matching
internal IDs), and "ambiguous" (other candidates with matching internal IDs). If "external" IDs were given as input, features based on them are initially classed as "unknown". Also in this case, a
support vector machine (SVM) is trained on the "positive" and "negative" candidates, to distinguish between the two classes based on the different OpenSWATH quality scores (plus an RT deviation
score). After parameter optimization by cross-validation, the resulting SVM is used to predict the probability of "unknown" feature candidates being positives.

<B>5. Feature filtering</B>

Feature candidates are filtered so that at most one feature per peptide and charge state remains.
For assays with internal IDs, only candidates previously classed as "positive" are kept.
For assays based solely on external IDs, the feature candidate with the highest SVM probability is selected and kept (possibly subject to the @p svm:min_prob threshold).

<B>6. Elution model fitting</B>

Elution models can be fitted to the features to improve the quantification.
For robustness, one model is fitted to all isotopic mass traces of a feature in parallel.
A symmetric (Gaussian) and an asymmetric (exponential-Gaussian hybrid) model type are available.
The fitted models are checked for plausibility before they are accepted.

Finally the results (feature maps, parameter @p out) are returned.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FeatureFinderIdentification.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FeatureFinderIdentification.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPFeatureFinderIdentification : public TOPPBase
{
public:
  // TODO
  // cppcheck-suppress uninitMemberVar
  TOPPFeatureFinderIdentification() :
      TOPPBase("FeatureFinderIdentification", "Detects features in MS1 data based on peptide identifications.", true,
               {{"Weisser H, Choudhary JS", "Targeted Feature Detection for Data-Dependent Shotgun Proteomics", "J. Proteome Res. 2017; 16, 8:2964-2974", "10.1021/acs.jproteome.7b00248"}})
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file: LC-MS raw data");
    setValidFormats_("in", {"mzML"});
    registerInputFile_("id", "<file>", "", "Input file: Peptide identifications derived directly from 'in'");
    setValidFormats_("id", {"idXML"});
    registerInputFile_("id_ext", "<file>", "", "Input file: 'External' peptide identifications (e.g. from aligned runs)", false);
    setValidFormats_("id_ext", {"idXML"});
    registerOutputFile_("out", "<file>", "", "Output file: Features");
    setValidFormats_("out", {"featureXML"});
    registerOutputFile_("lib_out", "<file>", "", "Output file: Assay library", false);
    setValidFormats_("lib_out", {"traML"});
    registerOutputFile_("chrom_out", "<file>", "", "Output file: Chromatograms", false);
    setValidFormats_("chrom_out", {"mzML"});
    registerOutputFile_("candidates_out", "<file>", "", "Output file: Feature candidates (before filtering and model fitting)", false);
    setValidFormats_("candidates_out", {"featureXML"});
    registerInputFile_("candidates_in", "<file>", "",
                       "Input file: Feature candidates from a previous run. If set, only feature classification and elution model fitting are carried out, if enabled. Many parameters are ignored.",
                       false, true);
    setValidFormats_("candidates_in", {"featureXML"});

    Param algo_with_subsection;
    Param subsection = FeatureFinderIdentificationAlgorithm().getDefaults();
    subsection.remove("candidates_out");
    algo_with_subsection.insert("", subsection);
    registerFullParam_(algo_with_subsection);
  }

  ExitCodes main_(int, const char**) override
  {
    FeatureMap features;

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String out = getStringOption_("out");
    String candidates_out = getStringOption_("candidates_out");

    String candidates_in = getStringOption_("candidates_in");

    FeatureFinderIdentificationAlgorithm ffid_algo;
    ffid_algo.getProgressLogger().setLogType(log_type_);
    ffid_algo.setParameters(getParam_().copySubset(FeatureFinderIdentificationAlgorithm().getDefaults()));

    if (candidates_in.empty())
    {
      String in = getStringOption_("in");
      String id = getStringOption_("id");
      String id_ext = getStringOption_("id_ext");
      String lib_out = getStringOption_("lib_out");
      String chrom_out = getStringOption_("chrom_out");

      //-------------------------------------------------------------
      // load input
      //-------------------------------------------------------------
      OPENMS_LOG_INFO << "Loading input data..." << endl;
      FileHandler mzml;
      mzml.getOptions().addMSLevel(1);
      mzml.loadExperiment(in, ffid_algo.getMSData(), {FileTypes::MZML}, log_type_);

      vector<PeptideIdentification> peptides, peptides_ext;
      vector<ProteinIdentification> proteins, proteins_ext;

      // "internal" IDs:
      FileHandler().loadIdentifications(id, proteins, peptides, {FileTypes::IDXML});

      // "external" IDs:
      if (!id_ext.empty())
      {
        FileHandler().loadIdentifications(id_ext, proteins_ext, peptides_ext, {FileTypes::IDXML});
      }

      //-------------------------------------------------------------
      // run feature detection
      //-------------------------------------------------------------
      ffid_algo.run(peptides, proteins, peptides_ext, proteins_ext, features, FeatureMap(), in);

      // write auxiliary output:
      bool keep_chromatograms = !chrom_out.empty();
      bool keep_library = !lib_out.empty();

      // keep assay data for output?
      if (keep_library)
      {
        FileHandler().storeTransitions(lib_out, ffid_algo.getLibrary(), {FileTypes::TRAML});
      }

      // keep chromatogram data for output?
      if (keep_chromatograms)
      {
        addDataProcessing_(ffid_algo.getChromatograms(), getProcessingInfo_(DataProcessing::FILTERING));
        FileHandler().storeExperiment(chrom_out, ffid_algo.getChromatograms(), {FileTypes::MZML});
        ffid_algo.getChromatograms().clear(true);
      }

      addDataProcessing_(features, getProcessingInfo_(DataProcessing::QUANTITATION));
    }
    else
    {
      //-------------------------------------------------------------
      // load feature candidates
      //-------------------------------------------------------------
      OPENMS_LOG_INFO << "Reading feature candidates from a previous run..." << endl;
      FileHandler().loadFeatures(candidates_in, features, {FileTypes::FEATUREXML});
      OPENMS_LOG_INFO << "Found " << features.size() << " feature candidates in total." << endl;
      ffid_algo.runOnCandidates(features);
    }

    //-------------------------------------------------------------
    // write output
    //-------------------------------------------------------------

    OPENMS_LOG_INFO << "Writing final results..." << endl;
    FileHandler().storeFeatures(out, features, {FileTypes::FEATUREXML});


    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
