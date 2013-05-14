// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Steffen Sass, Holger Plattfaut, Bastian Blank $
// --------------------------------------------------------------------------

//OpenMS includes
#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>
#include <OpenMS/METADATA/MSQuantifications.h>

#include <OpenMS/FILTERING/DATAREDUCTION/SILACAnalyzer.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <OpenMS/COMPARISON/CLUSTERING/SILACClustering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

//Contrib includes
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>
#include <locale>
#include <iomanip>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_SILACAnalyzer SILACAnalyzer

  @brief Identifies peptide pairs in LC-MS data and determines their relative abundance.

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ SILACAnalyzer \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDMapper</td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
    </tr>
  </table>
</CENTER>

  SILACAnalyzer is a tool for the fully automated analysis of quantitative proteomics data. It identifies pairs of isotopic envelopes with fixed m/z separation. It requires no prior sequence identification of the peptides. In what follows we first explain the algorithm and then discuss the tuning of its parameters.

  <b>Algorithm</b>

  The algorithm is divided into three parts: filtering, clustering and linear fitting, see Fig. (d), (e) and (f). In the following discussion let us consider a particular mass spectrum at retention time 1350 s, see Fig. (a). It contains a peptide of mass 1492 Da and its 6 Da heavier labelled counterpart. Both are doubly charged in this instance. Their isotopic envelopes therefore appear at 746 and 749 in the spectrum. The isotopic peaks within each envelope are separated by 0.5. The spectrum was recorded at finite intervals. In order to read accurate intensities at arbitrary m/z we spline-fit over the data, see Fig. (b).

  We would like to search for such peptide pairs in our LC-MS data set. As a warm-up let us consider a standard intensity cut-off filter, see Fig. (c). Scanning through the entire m/z range (red dot) only data points with intensities above a certain threshold pass the filter. Unlike such a local filter, the filter used in our algorithm takes intensities at a range of m/z positions into account, see Fig. (d). A data point (red dot) passes if
  - all six intensities at m/z, m/z+0.5, m/z+1, m/z+3, m/z+3.5 and m/z+4 lie above a certain threshold,
  - the intensity profiles in neighbourhoods around all six m/z positions show a good correlation and
  - the relative intensity ratios within a peptide agree up to a factor with the ratios of a theoretic averagine model.

  Let us now filter not only a single spectrum but all spectra in our data set. Data points that pass the filter form clusters in the t-m/z plane, see Fig. (e). Each cluster corresponds to the mono-isotopic mass trace of the lightest peptide of a SILAC pattern. We now use hierarchical clustering methods to assign each data point to a specific cluster. The optimum number of clusters is determined by maximizing the silhouette width of the partitioning. Each data point in a cluster corresponds to three pairs of intensities (at [m/z, m/z+3], [m/z+0.5, m/z+3.5] and [m/z+1, m/z+4]). A plot of all intensity pairs in a cluster shows a clear linear correlation, see Fig. (f). Using linear regression we can determine the relative amounts of labelled and unlabelled peptides in the sample.

  @image html SILACAnalyzer_algorithm.png

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_SILACAnalyzer.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_SILACAnalyzer.html

  <b>Parameter Tuning</b>

  SILACAnalyzer can detect SILAC patterns of any number of peptides, i.e. doublets (pairs), triplets, quadruplets et cetera.

  <i>input:</i>
  - in [*.mzML] - LC-MS dataset to be analyzed
  - ini [*.ini] - file containing all parameters (see discussion below)

  <i>standard output:</i>
  - out [*.consensusXML] - contains the list of identified peptides (retention time and m/z of the lightest peptide, ratios)

  <i>optional output:</i>
  - out_clusters [*.consensusXML] - contains the complete set of data points passing the filters, see Fig. (e)

  The results of an analysis can easily visualized within TOPPView. Simply load *.consensusXML and *.featureXML as layers over the original *.mzML.

  Parameters in section <i>algorithm:</i>
  - <i>allow_missing_peaks</i> - Low intensity peaks might be missing from the isotopic pattern of some of the peptides. Specify if such peptides should be included in the analysis.
  - <i>rt_threshold</i> - Upper bound for the retention time [s] over which a characteristic peptide elutes.
  - <i>rt_min</i> - Lower bound for the retentions time [s].
  - <i>intensity_cutoff</i> - Lower bound for the intensity of isotopic peaks in a SILAC pattern.
  - <i>intensity_correlation</i> - Lower bound for the Pearson correlation coefficient, which measures how well intensity profiles of different isotopic peaks correlate.
  - <i>model_deviation</i> - Upper bound on the factor by which the ratios of observed isotopic peaks are allowed to differ from the ratios of the theoretic averagine model, i.e. ( theoretic_ratio / model_deviation ) < observed_ratio < ( theoretic_ratio * model_deviation ).

  Parameters in section <i>sample:</i>
  - <i>labels</i> - Labels used for labelling the sample. [...] specifies the labels for a single sample. For example, [Lys4,Arg6][Lys8,Arg10] describes a mixtures of three samples. One of them unlabelled, one labelled with Lys4 and Arg6 and a third one with Lys8 and Arg10. For permitted labels see section <i>labels</i>.
  - <i>charge</i> - Range of charge states in the sample, i.e. min charge : max charge.
  - <i>missed_cleavages</i> - Maximum number of missed cleavages.
  - <i>peaks_per_peptide</i> - Range of peaks per peptide in the sample, i.e. min peaks per peptide : max peaks per peptide.

 Parameters in section <i>labels:</i>
 This section contains a list of all isotopic labels currently available for analysis of SILAC data with SILACAnalyzer.

 <b>References:</b>
  @n L. Nilse, M. Sturm, D. Trudgian, M. Salek, P. Sims, K. Carroll, S. Hubbard,  <a href="http://www.springerlink.com/content/u40057754100v71t">SILACAnalyzer - a tool for differential quantitation of stable isotope derived data</a>, in F. Masulli, L. Peterson, and R. Tagliaferri (Eds.): CIBB 2009, LNBI 6160, pp. 4555, 2010.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSILACAnalyzer :
  public TOPPBase
{
private:

  // input and output files
  String in;
  String out;
  String out_clusters;
  String out_features;
  String out_mzq;

  String out_filters;
  String in_filters;
  String out_debug;

  // section "sample"
  String selected_labels;
  UInt charge_min;
  UInt charge_max;
  Int missed_cleavages;
  UInt isotopes_per_peptide_min;
  UInt isotopes_per_peptide_max;

  // section "algorithm"
  DoubleReal rt_threshold;
  DoubleReal rt_min;
  DoubleReal intensity_cutoff;
  DoubleReal intensity_correlation;
  DoubleReal model_deviation;
  bool allow_missing_peaks;

  typedef SILACClustering Clustering;

public:
  TOPPSILACAnalyzer() :
    TOPPBase("SILACAnalyzer", "Determination of peak ratios in LC-MS data", true), allow_missing_peaks(true)
  {
  }

  //--------------------------------------------------
  // set structure of ini file
  //--------------------------------------------------

  void registerOptionsAndFlags_()
  {
    // create flag for input file (.mzML)
    registerInputFile_("in", "<file>", "", "Raw LC-MS data to be analyzed. (Profile data required. Will not work with centroided data!)");
    setValidFormats_("in", StringList::create("mzML"));
    // create flag for output file (.consensusXML)
    registerOutputFile_("out", "<file>", "", "Set of all identified peptide groups (i.e. peptide pairs or triplets or singlets or ..). The m/z-RT positions correspond to the lightest peptide in each group.", false);
    setValidFormats_("out", StringList::create("consensusXML"));
    // create optional flag for additional clusters output file (.featureXML)
    registerOutputFile_("out_clusters", "<file>", "", "Optional debug output containing data points passing all filters, hence belonging to a SILAC pattern. Points of the same colour correspond to the mono-isotopic peak of the lightest peptide in a pattern.", false, true);
    setValidFormats_("out_clusters", StringList::create("consensusXML"));
    registerOutputFile_("out_features", "<file>", "", "Optional output file containing the individual peptide features in \'out\'.", false, true);
    setValidFormats_("out_features", StringList::create("featureXML"));
    registerOutputFile_("out_mzq", "<file>", "", "Optional output file of MzQuantML.", false, true);
    setValidFormats_("out_mzq", StringList::create("mzq"));

    // create optional flag for additional output file (.consensusXML) to store filter results
    registerOutputFile_("out_filters", "<file>", "", "Optional output file containing all points that passed the filters as txt. Suitable as input for \'in_filters\' to perform clustering without preceding filtering process.", false, true);
    setValidFormats_("out_filters", StringList::create("consensusXML"));
    // create optional flag for additional input file (.consensusXML) to load filter results
    registerInputFile_("in_filters", "<file>", "", "Optional input file containing all points that passed the filters as txt. Use output from \'out_filters\' to perform clustering only.", false, true);
    setValidFormats_("in_filters", StringList::create("consensusXML"));
    registerStringOption_("out_debug", "<filebase>", "", "Filename base for debug output.", false, true);

    // create section "labels" for adjusting masses of labels
    registerSubsection_("labels", "Isotopic labels that can be specified in section \'sample\'.");
    // create section "sample" for adjusting sample parameters
    registerSubsection_("sample", "Parameters describing the sample and its labels.");
    // create section "algorithm" for adjusting algorithm parameters
    registerSubsection_("algorithm", "Parameters for the algorithm.");

    // create flag for missing peaks
    registerFlag_("algorithm:allow_missing_peaks", "Low intensity peaks might be missing from the isotopic pattern of some of the peptides. Should such peptides be included in the analysis?", true);
  }

  // create prameters for sections (set default values and restrictions)
  Param getSubsectionDefaults_(const String & section) const
  {
    Param defaults;


    //--------------------------------------------------
    // section labels
    //--------------------------------------------------

    if (section == "labels")
    {
      // create labels that can be chosen in section "sample/labels"
      defaults.setValue("Arg6", 6.0201290268, "Arg6 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Arg6", 0.0);
      defaults.setValue("Arg10", 10.008268600, "Arg10 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Arg10", 0.0);
      defaults.setValue("Lys4", 4.0251069836, "Lys4 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Lys4", 0.0);
      defaults.setValue("Lys6", 6.0201290268, "Lys6 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Lys6", 0.0);
      defaults.setValue("Lys8", 8.0141988132, "Lys8 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Lys8", 0.0);
      defaults.setValue("dICPL4", 4.025107, "mass difference between isotope-coded protein labels ICPL 4 and ICPL 0", StringList::create("advanced"));
      defaults.setMinFloat("dICPL4", 0.0);
      defaults.setValue("dICPL6", 6.020129, "mass difference between isotope-coded protein labels ICPL 6 and ICPL 0", StringList::create("advanced"));
      defaults.setMinFloat("dICPL6", 0.0);
      defaults.setValue("dICPL10", 10.045236, "mass difference between isotope-coded protein labels ICPL 10 and ICPL 0", StringList::create("advanced"));
      defaults.setMinFloat("dICPL10", 0.0);
      defaults.setValue("Methyl4", 4.0202, "Methyl4 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Methyl4", 0.0);
      defaults.setValue("Methyl8", 8.0202, "Methyl8 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Methyl8", 0.0);
      defaults.setValue("Methyl12", 12.0202, "Methyl12 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Methyl12", 0.0);
      defaults.setValue("Methyl16", 16.0202, "Methyl16 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Methyl16", 0.0);
      defaults.setValue("Methyl24", 24.0202, "Methyl24 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Methyl24", 0.0);
      defaults.setValue("Methyl32", 32.0202, "Methyl32 mass shift", StringList::create("advanced"));
      defaults.setMinFloat("Methyl32", 0.0);
    }


    //--------------------------------------------------
    // section sample
    //--------------------------------------------------

    if (section == "sample")
    {
      defaults.setValue("labels", "[Lys8,Arg10]", "Labels used for labelling the sample. [...] specifies the labels for a single sample. For example, [Lys4,Arg6][Lys8,Arg10] describes a mixtures of three samples. One of them unlabelled, one labelled with Lys4 and Arg6 and a third one with Lys8 and Arg10. For permitted labels see \'advanced parameters\', section \'labels\'. If left empty the tool identifies singlets, i.e. acts as peptide feature finder.");
      defaults.setValue("charge", "2:4", "Range of charge states in the sample, i.e. min charge : max charge.");
      defaults.setValue("missed_cleavages", 0, "Maximum number of missed cleavages.");
      defaults.setMinInt("missed_cleavages", 0);
      defaults.setValue("peaks_per_peptide", "3:5", "Range of peaks per peptide in the sample, i.e. min peaks per peptide : max peaks per peptide. For example 3:6, if isotopic peptide patterns in the sample consist of either three, four, five or six isotopic peaks. ", StringList::create("advanced"));
    }


    //--------------------------------------------------
    // section algorithm
    //--------------------------------------------------

    if (section == "algorithm")
    {
      defaults.setValue("rt_threshold", 30.0, "Typical retention time [s] over which a characteristic peptide elutes. (This is not an upper bound. Peptides that elute for longer will be reported.)");
      defaults.setMinFloat("rt_threshold", 0.0);
      defaults.setValue("rt_min", 0.0, "Lower bound for the retention time [s].");
      defaults.setMinFloat("rt_min", 0.0);
      defaults.setValue("intensity_cutoff", 1000.0, "Lower bound for the intensity of isotopic peaks in a SILAC pattern.");
      defaults.setMinFloat("intensity_cutoff", 0.0);
      defaults.setValue("intensity_correlation", 0.7, "Lower bound for the Pearson correlation coefficient, which measures how well intensity profiles of different isotopic peaks correlate.");
      defaults.setMinFloat("intensity_correlation", 0.0);
      defaults.setMaxFloat("intensity_correlation", 1.0);
      defaults.setValue("model_deviation", 3.0, "Upper bound on the factor by which the ratios of observed isotopic peaks are allowed to differ from the ratios of the theoretic averagine model, i.e. ( theoretic_ratio / model_deviation ) < observed_ratio < ( theoretic_ratio * model_deviation ).");
      defaults.setMinFloat("model_deviation", 1.0);
    }

    return defaults;
  }

  //--------------------------------------------------
  // handle parameters (read in and format given parameters)
  //--------------------------------------------------

  void handleParameters_sample()
  {

    //--------------------------------------------------
    // section sample
    //--------------------------------------------------

    // get selected labels
    selected_labels = getParam_().getValue("sample:labels");

    // get selected missed_cleavages
    missed_cleavages = getParam_().getValue("sample:missed_cleavages");

    // get selected charge range
    String charge_string = getParam_().getValue("sample:charge");
    DoubleReal charge_min_temp, charge_max_temp;
    parseRange_(charge_string, charge_min_temp, charge_max_temp);
    charge_min = charge_min_temp;
    charge_max = charge_max_temp;

    // check if charge_min is smaller than charge max, if not swap
    if (charge_min > charge_max)
      swap(charge_min, charge_max);

    // get selected peaks range
    String isotopes_per_peptide_string = getParam_().getValue("sample:peaks_per_peptide");
    DoubleReal isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp;
    parseRange_(isotopes_per_peptide_string, isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp);
    isotopes_per_peptide_min = isotopes_per_peptide_min_temp;
    isotopes_per_peptide_max = isotopes_per_peptide_max_temp;

    //check if isotopes_per_peptide_min is smaller than isotopes_per_peptide_max, if not swap
    if (isotopes_per_peptide_min > isotopes_per_peptide_max)
      swap(isotopes_per_peptide_min, isotopes_per_peptide_max);
  }

  void handleParameters_algorithm()
  {
    //--------------------------------------------------
    // section algorithm
    //--------------------------------------------------

    rt_threshold = getParam_().getValue("algorithm:rt_threshold");
    rt_min = getParam_().getValue("algorithm:rt_min");
    intensity_cutoff = getParam_().getValue("algorithm:intensity_cutoff");
    intensity_correlation = getParam_().getValue("algorithm:intensity_correlation");
    model_deviation = getParam_().getValue("algorithm:model_deviation");
    allow_missing_peaks = getFlag_("algorithm:allow_missing_peaks");
  }

  void handleParameters_labels(map<String, DoubleReal> & label_identifiers)
  {

    //--------------------------------------------------
    // section labels
    //--------------------------------------------------

    // create map of pairs (label as string, mass shift as double)
    label_identifiers.insert(make_pair("Arg6", getParam_().getValue("labels:Arg6")));
    label_identifiers.insert(make_pair("Arg10", getParam_().getValue("labels:Arg10")));
    label_identifiers.insert(make_pair("Lys4", getParam_().getValue("labels:Lys4")));
    label_identifiers.insert(make_pair("Lys6", getParam_().getValue("labels:Lys6")));
    label_identifiers.insert(make_pair("Lys8", getParam_().getValue("labels:Lys8")));
    label_identifiers.insert(make_pair("Methyl4", getParam_().getValue("labels:Methyl4")));
    label_identifiers.insert(make_pair("Methyl8", getParam_().getValue("labels:Methyl8")));
    label_identifiers.insert(make_pair("Methyl12", getParam_().getValue("labels:Methyl12")));
    label_identifiers.insert(make_pair("Methyl16", getParam_().getValue("labels:Methyl16")));
    label_identifiers.insert(make_pair("Methyl24", getParam_().getValue("labels:Methyl24")));
    label_identifiers.insert(make_pair("Methyl32", getParam_().getValue("labels:Methyl32")));
    label_identifiers.insert(make_pair("dICPL4", getParam_().getValue("labels:dICPL4")));
    label_identifiers.insert(make_pair("dICPL6", getParam_().getValue("labels:dICPL6")));
    label_identifiers.insert(make_pair("dICPL10", getParam_().getValue("labels:dICPL10")));

  }

  void handleParameters()
  {
    // get input file (.mzML)
    in = getStringOption_("in");
    // get name of output file (.consensusXML)
    out = getStringOption_("out");
    // get name of additional clusters output file (.consensusXML)
    out_clusters = getStringOption_("out_clusters");
    out_features = getStringOption_("out_features");
    out_mzq = getStringOption_("out_mzq");

    // get name of additional filters output file (.consensusXML)
    out_filters = getStringOption_("out_filters");
    // get name of additional filters input file (.consensusXML)
    in_filters = getStringOption_("in_filters");
    out_debug = getStringOption_("out_debug");
  }

  //--------------------------------------------------
  // filtering
  //--------------------------------------------------

  ExitCodes main_(int, const char **)
  {
    // data to be passed through the algorithm
    vector<vector<SILACPattern> > data;
    MSQuantifications msq;
    vector<Clustering *> cluster_data;

    // 
    // Parameter handling
    // 
    map<String, DoubleReal> label_identifiers;   // list defining the mass shifts of each label (e.g. "Arg6" => 6.0201290268)
    handleParameters_sample();
    handleParameters_algorithm();
    handleParameters_labels(label_identifiers);
    handleParameters();

    // 
    // Initializing the SILACAnalzer with our parameters
    // 
    SILACAnalyzer analyzer;
    analyzer.setLogType(log_type_);
    analyzer.initialize(
      // section "sample"
      selected_labels,
      charge_min,
      charge_max,
      missed_cleavages,
      isotopes_per_peptide_min,
      isotopes_per_peptide_max,
      // section "algorithm"
      rt_threshold,
      rt_min,
      intensity_cutoff,
      intensity_correlation,
      model_deviation,
      allow_missing_peaks,
      // labels
      label_identifiers);


    //--------------------------------------------------
    // loading input from .mzML
    //--------------------------------------------------

    MzMLFile file;
    MSExperiment<Peak1D> exp;

    // only read MS1 spectra ...
    /*
    std::vector<int> levels;
    levels.push_back(1);
    file.getOptions().setMSLevels(levels);
    */
    file.setLogType(log_type_);
    file.load(in, exp);

    // set size of input map
    exp.updateRanges();

    // extract level 1 spectra
    exp.erase(remove_if(exp.begin(), exp.end(), InMSLevelRange<MSExperiment<Peak1D>::SpectrumType>(IntList::create("1"), true)), exp.end());

    // sort according to RT and MZ
    exp.sortSpectra();

    if (out_mzq != "")
    {
      vector<vector<String> > SILAClabels = analyzer.getSILAClabels(); // list of SILAC labels, e.g. selected_labels="[Lys4,Arg6][Lys8,Arg10]" => SILAClabels[0][1]="Arg6"

      std::vector<std::vector<std::pair<String, DoubleReal> > > labels;
      //add none label
      labels.push_back(std::vector<std::pair<String, DoubleReal> >(1, std::make_pair<String, DoubleReal>(String("none"), DoubleReal(0))));
      for (Size i = 0; i < SILAClabels.size(); ++i)       //SILACLabels MUST be in weight order!!!
      {
        std::vector<std::pair<String, DoubleReal> > one_label;
        for (UInt j = 0; j < SILAClabels[i].size(); ++j)
        {
          one_label.push_back(*(label_identifiers.find(SILAClabels[i][j])));              // this dereferencing would break if all SILAClabels would not have been checked before!
        }
        labels.push_back(one_label);
      }
      msq.registerExperiment(exp, labels);       //add assays
      msq.assignUIDs();
    }
    MSQuantifications::QUANT_TYPES quant_type = MSQuantifications::MS1LABEL;
    msq.setAnalysisSummaryQuantType(quant_type);    //add analysis_summary_

    //--------------------------------------------------
    // estimate peak width
    //--------------------------------------------------

    PeakWidthEstimator::Result peak_width;
    try
    {
      peak_width = analyzer.estimatePeakWidth(exp);
    }
    catch (Exception::InvalidSize &)
    {
      writeLog_("Error: Unable to estimate peak width of input data.");
      return INCOMPATIBLE_INPUT_DATA;
    }


    if (in_filters == "")
    {
      //--------------------------------------------------
      // filter input data
      //--------------------------------------------------

      analyzer.filterData(exp, peak_width, data); 

      //--------------------------------------------------
      // store filter results
      //--------------------------------------------------

      if (out_filters != "")
      {
        ConsensusMap map;
        for (std::vector<std::vector<SILACPattern> >::const_iterator it = data.begin(); it != data.end(); ++it)
        {
          analyzer.generateFilterConsensusByPattern(map, *it);
        }
        analyzer.writeConsensus(out_filters, map);
      }
    }
    else
    {
      //--------------------------------------------------
      // load filter results
      //--------------------------------------------------

      ConsensusMap map;
      analyzer.readConsensus(in_filters, map);
      analyzer.readFilterConsensusByPattern(map, data);
    }

    //--------------------------------------------------
    // clustering
    //--------------------------------------------------

    analyzer.clusterData(exp, peak_width, cluster_data, data);

    //--------------------------------------------------------------
    // write output
    //--------------------------------------------------------------

    if (out_debug != "")
    {
      std::ofstream out((out_debug + ".clusters.csv").c_str());

      vector<vector<DoubleReal> > massShifts = analyzer.getMassShifts(); // list of mass shifts

      // generate header
      out
      << std::fixed << std::setprecision(8)
      << "ID,RT,MZ_PEAK,CHARGE";
      for (UInt i = 1; i <= massShifts[0].size(); ++i)
      {
        out << ",DELTA_MASS_" << i + 1;
      }
      for (UInt i = 0; i <= massShifts[0].size(); ++i)
      {
        for (UInt j = 1; j <= isotopes_per_peptide_max; ++j)
        {
          out << ",INT_PEAK_" << i + 1 << '_' << j;
        }
      }
      out << ",MZ_RAW";
      for (UInt i = 0; i <= massShifts[0].size(); ++i)
      {
        for (UInt j = 1; j <= isotopes_per_peptide_max; ++j)
        {
          out << ",INT_RAW_" << i + 1 << '_' << j;
        }
      }
      for (UInt i = 0; i <= massShifts[0].size(); ++i)
      {
        for (UInt j = 1; j <= isotopes_per_peptide_max; ++j)
        {
          out << ",MZ_RAW_" << i + 1 << '_' << j;
        }
      }
      out << '\n';

      // write data
      UInt cluster_id = 0;
      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        analyzer.generateClusterDebug(out, **it, cluster_id);
      }
    }

    if (out != "")
    {
      ConsensusMap map;

      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        analyzer.generateClusterConsensusByCluster(map, **it);
      }

      // XXX: Need a map per mass shift
      ConsensusMap::FileDescriptions & desc = map.getFileDescriptions();
      UInt id = 0;
      for (ConsensusMap::FileDescriptions::iterator it = desc.begin(); it != desc.end(); ++it, ++id)
      {
        if (!test_mode_) it->second.filename = in;
        // XXX: Write correct label
        // it->second.label = id;
      }


      std::set<DataProcessing::ProcessingAction> actions;
      actions.insert(DataProcessing::DATA_PROCESSING);
      actions.insert(DataProcessing::PEAK_PICKING);
      actions.insert(DataProcessing::FILTERING);
      actions.insert(DataProcessing::QUANTITATION);

      addDataProcessing_(map, getProcessingInfo_(actions));


      analyzer.writeConsensus(out, map);
      if (out_mzq != "")
      {
        ConsensusMap numap(map);
        //calc. ratios
        for (ConsensusMap::iterator cit = numap.begin(); cit != numap.end(); ++cit)
        {
          //~ make ratio templates
          std::vector<ConsensusFeature::Ratio> rts;
          for (std::vector<MSQuantifications::Assay>::const_iterator ait = msq.getAssays().begin() + 1; ait != msq.getAssays().end(); ++ait)
          {
            ConsensusFeature::Ratio r;
            r.numerator_ref_ = String(msq.getAssays().begin()->uid_);
            r.denominator_ref_ = String(ait->uid_);
            r.description_.push_back("Simple ratio calc");
            r.description_.push_back("light to medium/.../heavy");
            //~ "<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001132\" name=\"peptide ratio\"/>"
            rts.push_back(r);
          }
          const std::set<FeatureHandle, FeatureHandle::IndexLess> & feature_handles = cit->getFeatures();
          if (feature_handles.size() > 1)
          {
            std::set<FeatureHandle, FeatureHandle::IndexLess>::const_iterator fit = feature_handles.begin();             // this is unlabeled
            fit++;
            for (; fit != feature_handles.end(); ++fit)
            {
              Size ri = std::distance(feature_handles.begin(), fit);
              rts[ri - 1].ratio_value_ =  feature_handles.begin()->getIntensity() / fit->getIntensity();             // a proper silacalanyzer algo should never have 0-intensities so no 0devison ...
            }
          }

          cit->setRatios(rts);
        }
        msq.addConsensusMap(numap);        //add SILACAnalyzer result

        //~ msq.addFeatureMap();//add SILACAnalyzer evidencetrail as soon as clear what is realy contained in the featuremap
        //~ add AuditCollection - no such concept in TOPPTools yet
        analyzer.writeMzQuantML(out_mzq, msq);
      }
    }

    if (out_clusters != "")
    {
      ConsensusMap map;
      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        UInt cluster_id = 0;
        analyzer.generateClusterConsensusByPattern(map, **it, cluster_id);
      }

      ConsensusMap::FileDescription & desc = map.getFileDescriptions()[0];
      desc.filename = in;
      desc.label = "Cluster";

      analyzer.writeConsensus(out_clusters, map);
    }

    if (out_features != "")
    {
      FeatureMap<> map;
      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        analyzer.generateClusterFeatureByCluster(map, **it);
      }

      analyzer.writeFeatures(out_features, map);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPSILACAnalyzer tool;
  return tool.main(argc, argv);
}

//@endcond

