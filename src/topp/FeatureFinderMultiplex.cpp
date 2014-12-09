// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Authors: Lars Nilse, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>

#include <OpenMS/METADATA/MSQuantifications.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/COMPARISON/CLUSTERING/GridBasedCluster.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

//Contrib includes
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <QDir>

//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>
#include <locale>
#include <iomanip>

using namespace std;
using namespace OpenMS;
using namespace boost::math;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_FeatureFinderMultiplex FeatureFinderMultiplex

  @brief Identifies peptide pairs in LC-MS data and determines their relative abundance.

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderMultiplex \f$ \longrightarrow \f$</td>
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

  FeatureFinderMultiplex is a tool for the fully automated analysis of quantitative proteomics data. It identifies pairs of isotopic envelopes with fixed m/z separation. It requires no prior sequence identification of the peptides. In what follows we first explain the algorithm and then discuss the tuning of its parameters.

  <b>Algorithm</b>

  The algorithm is divided into three parts: filtering, clustering and linear fitting, see Fig. (d), (e) and (f). In the following discussion let us consider a particular mass spectrum at retention time 1350 s, see Fig. (a). It contains a peptide of mass 1492 Da and its 6 Da heavier labelled counterpart. Both are doubly charged in this instance. Their isotopic envelopes therefore appear at 746 and 749 in the spectrum. The isotopic peaks within each envelope are separated by 0.5. The spectrum was recorded at finite intervals. In order to read accurate intensities at arbitrary m/z we spline-fit over the data, see Fig. (b).

  We would like to search for such peptide pairs in our LC-MS data set. As a warm-up let us consider a standard intensity cut-off filter, see Fig. (c). Scanning through the entire m/z range (red dot) only data points with intensities above a certain threshold pass the filter. Unlike such a local filter, the filter used in our algorithm takes intensities at a range of m/z positions into account, see Fig. (d). A data point (red dot) passes if
  - all six intensities at m/z, m/z+0.5, m/z+1, m/z+3, m/z+3.5 and m/z+4 lie above a certain threshold,
  - the intensity profiles in neighbourhoods around all six m/z positions show a good correlation and
  - the relative intensity ratios within a peptide agree up to a factor with the ratios of a theoretic averagine model.

  Let us now filter not only a single spectrum but all spectra in our data set. Data points that pass the filter form clusters in the t-m/z plane, see Fig. (e). Each cluster corresponds to the mono-isotopic mass trace of the lightest peptide of a SILAC pattern. We now use hierarchical clustering methods to assign each data point to a specific cluster. The optimum number of clusters is determined by maximizing the silhouette width of the partitioning. Each data point in a cluster corresponds to three pairs of intensities (at [m/z, m/z+3], [m/z+0.5, m/z+3.5] and [m/z+1, m/z+4]). A plot of all intensity pairs in a cluster shows a clear linear correlation, see Fig. (f). Using linear regression we can determine the relative amounts of labelled and unlabelled peptides in the sample.

  @image html SILACAnalyzer_algorithm.png

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_FeatureFinderMultiplex.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_FeatureFinderMultiplex.html

  <b>Parameter Tuning</b>

  FeatureFinderMultiplex can detect isotope patterns of any number of peptides, i.e. doublets (pairs), triplets, quadruplets et cetera.

  <i>input:</i>
  - in [*.mzML] - LC-MS dataset to be analyzed
  - ini [*.ini] - file containing all parameters (see discussion below)

  <i>output:</i>
  - out [*.consensusXML] - contains the list of identified peptide multiples (retention time and m/z of the lightest peptide, ratios)
  - out_features [*.featureXML] - contains the list of individual peptides
  - out_mzq [*.mzq] - contains the results in mzQuantML format

  The results of an analysis can easily visualized within TOPPView. Simply load *.consensusXML and *.featureXML as layers over the original *.mzML.

  Parameters in section <i>algorithm:</i>
  - <i>allow_missing_peaks</i> - Low intensity peaks might be missing from the isotopic pattern of some of the peptides. Specify if such peptides should be included in the analysis.
  - <i>rt_typical</i> - Upper bound for the retention time [s] over which a characteristic peptide elutes.
  - <i>rt_min</i> - Lower bound for the retentions time [s].
  - <i>intensity_cutoff</i> - Lower bound for the intensity of isotopic peaks in a peptide pattern.
  - <i>peptide_similarity</i> - Lower bound for the Pearson correlation coefficient, which measures how well intensity profiles of different isotopic peaks correlate.
  - <i>averagine_similarity</i> - Lower bound for the Pearson correlation coefficient, which measures how well the isotope patterns match the theoretical averagine model.

  Parameters in section <i>algorithm:</i>
  - <i>labels</i> - Labels used for labelling the sample. [...] specifies the labels for a single sample. For example, [][Lys4,Arg6][Lys8,Arg10] describes a mixtures of three samples. One of them unlabelled, one labelled with Lys4 and Arg6 and a third one with Lys8 and Arg10. For permitted labels see section <i>labels</i>.
  - <i>charge</i> - Range of charge states in the sample, i.e. min charge : max charge.
  - <i>missed_cleavages</i> - Maximum number of missed cleavages.
  - <i>isotopes_per_peptide</i> - Range of peaks per peptide in the sample, i.e. min peaks per peptide : max peaks per peptide.

 Parameters in section <i>labels:</i>
 This section contains a list of all isotopic labels currently available for analysis with FeatureFinderMultiplex.

 <b>References:</b>
  @n L. Nilse, M. Sturm, D. Trudgian, M. Salek, P. Sims, K. Carroll, S. Hubbard,  <a href="http://www.springerlink.com/content/u40057754100v71t">SILACAnalyzer - a tool for differential quantitation of stable isotope derived data</a>, in F. Masulli, L. Peterson, and R. Tagliaferri (Eds.): CIBB 2009, LNBI 6160, pp. 4555, 2010.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMultiplex :
  public TOPPBase
{
private:

  // input and output files
  String in_;
  String out_;
  String out_features_;
  String out_mzq_;
  String out_debug_;
  bool debug_;

  // section "algorithm"
  String labels_;
  std::vector<std::vector<String> > samples_labels_;
  unsigned charge_min_;
  unsigned charge_max_;
  unsigned missed_cleavages_;
  unsigned isotopes_per_peptide_min_;
  unsigned isotopes_per_peptide_max_;
  double rt_typical_;
  double rt_min_;
  double mz_tolerance_;
  bool mz_unit_; // ppm (true), Da (false)
  double intensity_cutoff_;
  double peptide_similarity_;
  double averagine_similarity_;
  double averagine_similarity_scaling_;
  bool knock_out_;

  // section "labels"
  map<String, double> label_massshift_;

public:
  TOPPFeatureFinderMultiplex() :
    TOPPBase("FeatureFinderMultiplex", "Determination of peak ratios in LC-MS data", true),
    debug_(false), charge_min_(1), charge_max_(1), missed_cleavages_(0), isotopes_per_peptide_min_(1), isotopes_per_peptide_max_(1), rt_typical_(0.0), rt_min_(0.0),
    mz_tolerance_(0.0), mz_unit_(true), intensity_cutoff_(0.0), peptide_similarity_(0.0), averagine_similarity_(0.0), averagine_similarity_scaling_(0.0), knock_out_(false)
  {
  }

  typedef std::vector<double> MassPattern; // list of mass shifts

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Raw LC-MS data to be analyzed. (Profile data required. Will not work with centroided data!)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Set of all identified peptide groups (i.e. peptide pairs or triplets or singlets or ..). The m/z-RT positions correspond to the lightest peptide in each group.", false);
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out_features", "<file>", "", "Optional output file containing the individual peptide features in \'out\'.", false, true);
    setValidFormats_("out_features", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out_mzq", "<file>", "", "Optional output file of MzQuantML.", false, true);
    setValidFormats_("out_mzq", ListUtils::create<String>("mzq"));
    registerStringOption_("out_debug", "<out_dir>", "", "Directory for debug output.", false, true);

    registerSubsection_("algorithm", "Parameters for the algorithm.");
    registerSubsection_("labels", "Isotopic labels that can be specified in section \'algorithm:labels\'.");

  }

  // create parameters for sections (set default values and restrictions)
  Param getSubsectionDefaults_(const String& section) const
  {
    Param defaults;

    if (section == "algorithm")
    {
      defaults.setValue("labels", "[][Lys8,Arg10]", "Labels used for labelling the samples. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL");
      defaults.setValue("charge", "1:4", "Range of charge states in the sample, i.e. min charge : max charge.");
      defaults.setValue("isotopes_per_peptide", "3:6", "Range of isotopes per peptide in the sample. For example 3:6, if isotopic peptide patterns in the sample consist of either three, four, five or six isotopic peaks. ", ListUtils::create<String>("advanced"));
      defaults.setValue("rt_typical", 40.0, "Typical retention time [s] over which a characteristic peptide elutes. (This is not an upper bound. Peptides that elute for longer will be reported.)");
      defaults.setMinFloat("rt_typical", 0.0);
      defaults.setValue("rt_min", 2.0, "Lower bound for the retention time [s]. (Any peptides seen for a shorter time period are not reported.)");
      defaults.setMinFloat("rt_min", 0.0);
      defaults.setValue("mz_tolerance", 6.0, "m/z tolerance for search of peak patterns.");
      defaults.setMinFloat("mz_tolerance", 0.0);
      defaults.setValue("mz_unit", "ppm", "Unit of the 'mz_tolerance' parameter.");
      defaults.setValidStrings("mz_unit", ListUtils::create<String>("Da,ppm"));
      defaults.setValue("intensity_cutoff", 1000.0, "Lower bound for the intensity of isotopic peaks.");
      defaults.setMinFloat("intensity_cutoff", 0.0);
      defaults.setValue("peptide_similarity", 0.5, "Two peptides in a multiplet are expected to have the same isotopic pattern. This parameter is a lower bound on their similarity.");
      defaults.setMinFloat("peptide_similarity", 0.0);
      defaults.setMaxFloat("peptide_similarity", 1.0);
      defaults.setValue("averagine_similarity", 0.4, "The isotopic pattern of a peptide should resemble the averagine model at this m/z position. This parameter is a lower bound on similarity between measured isotopic pattern and the averagine model.");
      defaults.setMinFloat("averagine_similarity", 0.0);
      defaults.setMaxFloat("averagine_similarity", 1.0);
      defaults.setValue("averagine_similarity_scaling", 0.75, "Let x denote this scaling factor, and p the averagine similarity parameter. For the detection of single peptides, the averagine parameter p is replaced by p' = p + x(1-p), i.e. x = 0 -> p' = p and x = 1 -> p' = 1. (For knock_out = true, peptide doublets and singlets are detected simulataneously. For singlets, the peptide similarity filter is irreleavant. In order to compensate for this 'missing filter', the averagine parameter p is replaced by the more restrictive p' when searching for singlets.)", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("averagine_similarity_scaling", 0.0);
      defaults.setMaxFloat("averagine_similarity_scaling", 1.0);
      defaults.setValue("missed_cleavages", 0, "Maximum number of missed cleavages due to incomplete digestion.");
      defaults.setMinInt("missed_cleavages", 0);
      defaults.setValue("knock_out", "false", "Is it likely that knock-outs are present? (Supported for doublex, triplex and quadruplex experiments only.)", ListUtils::create<String>("advanced"));
      defaults.setValidStrings("knock_out", ListUtils::create<String>("true,false"));
    }

    if (section == "labels")
    {
      // create labels that can be chosen in section "algorithm/labels"
      defaults.setValue("Arg6", 6.0201290268, "Label:13C(6)  |  C(-6) 13C(6)  |  unimod #188", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Arg6", 0.0);
      defaults.setValue("Arg10", 10.008268600, "Label:13C(6)15N(4)  |  C(-6) 13C(6) N(-4) 15N(4)  |  unimod #267", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Arg10", 0.0);
      defaults.setValue("Lys4", 4.0251069836, "Label:2H(4)  |  H(-4) 2H(4)  |  unimod #481", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Lys4", 0.0);
      defaults.setValue("Lys6", 6.0201290268, "Label:13C(6)  |  C(-6) 13C(6)  |  unimod #188", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Lys6", 0.0);
      defaults.setValue("Lys8", 8.0141988132, "Label:13C(6)15N(2)  |  C(-6) 13C(6) N(-2) 15N(2)  |  unimod #259", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Lys8", 0.0);
      defaults.setValue("Dimethyl0", 28.031300, "Dimethyl  |  H(4) C(2)  |  unimod #36", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Dimethyl0", 0.0);
      defaults.setValue("Dimethyl4", 32.056407, "Dimethyl:2H(4)  |  2H(4) C(2)  |  unimod #199", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Dimethyl4", 0.0);
      defaults.setValue("Dimethyl6", 34.063117, "Dimethyl:2H(4)13C(2)  |  2H(4) 13C(2)  |  unimod #510", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Dimethyl6", 0.0);
      defaults.setValue("Dimethyl8", 36.075670, "Dimethyl:2H(6)13C(2)  |  H(-2) 2H(6) 13C(2)  |  unimod #330", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("Dimethyl8", 0.0);
      defaults.setValue("ICPL0", 105.021464, "ICPL  |  H(3) C(6) N O  |  unimod #365", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("ICPL0", 0.0);
      defaults.setValue("ICPL4", 109.046571, "ICPL:2H(4)  |  H(-1) 2H(4) C(6) N O  |  unimod #687", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("ICPL4", 0.0);
      defaults.setValue("ICPL6", 111.041593, "ICPL:13C(6)  |  H(3) 13C(6) N O  |  unimod #364", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("ICPL6", 0.0);
      defaults.setValue("ICPL10", 115.066700, "ICPL:13C(6)2H(4)  |  H(-1) 2H(4) 13C(6) N O  |  unimod #866", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("ICPL10", 0.0);
      //defaults.setValue("18O", 2.004246, "Label:18O(1)  |  O(-1) 18O  |  unimod #258", ListUtils::create<String>("advanced"));
      //defaults.setMinFloat("18O", 0.0);
    }

    return defaults;
  }

  /**
   * @brief process parameters of 'input/output' section
   */
  void getParameters_in_out_()
  {
    in_ = getStringOption_("in");
    out_ = getStringOption_("out");
    out_features_ = getStringOption_("out_features");
    out_mzq_ = getStringOption_("out_mzq");
    out_debug_ = getStringOption_("out_debug");
    debug_ = !out_debug_.empty();
  }

  /**
   * @brief process parameters of 'algorithm' section
   */
  void getParameters_algorithm_()
  {
    // get selected labels
    labels_ = getParam_().getValue("algorithm:labels");
    samples_labels_ = splitLabelString_();

    // get selected charge range
    String charge_string = getParam_().getValue("algorithm:charge");
    double charge_min_temp, charge_max_temp;
    parseRange_(charge_string, charge_min_temp, charge_max_temp);
    charge_min_ = charge_min_temp;
    charge_max_ = charge_max_temp;
    if (charge_min_ > charge_max_)
    {
      swap(charge_min_, charge_max_);
    }

    // get isotopes per peptide range
    String isotopes_per_peptide_string = getParam_().getValue("algorithm:isotopes_per_peptide");
    double isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp;
    parseRange_(isotopes_per_peptide_string, isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp);
    isotopes_per_peptide_min_ = isotopes_per_peptide_min_temp;
    isotopes_per_peptide_max_ = isotopes_per_peptide_max_temp;
    if (isotopes_per_peptide_min_ > isotopes_per_peptide_max_)
    {
      swap(isotopes_per_peptide_min_, isotopes_per_peptide_max_);
    }

    rt_typical_ = getParam_().getValue("algorithm:rt_typical");
    rt_min_ = getParam_().getValue("algorithm:rt_min");
    mz_tolerance_ = getParam_().getValue("algorithm:mz_tolerance");
    mz_unit_ = (getParam_().getValue("algorithm:mz_unit") == "ppm");
    intensity_cutoff_ = getParam_().getValue("algorithm:intensity_cutoff");
    peptide_similarity_ = getParam_().getValue("algorithm:peptide_similarity");
    averagine_similarity_ = getParam_().getValue("algorithm:averagine_similarity");
    averagine_similarity_scaling_ = getParam_().getValue("algorithm:averagine_similarity_scaling");
    missed_cleavages_ = getParam_().getValue("algorithm:missed_cleavages");
    knock_out_ = (getParam_().getValue("algorithm:knock_out") == "true");
  }

  /**
   * @brief process parameters of 'labels' section
   */
  void getParameters_labels_()
  {
    // create map of pairs (label as string, mass shift as double)
    label_massshift_.insert(make_pair("Arg6", getParam_().getValue("labels:Arg6")));
    label_massshift_.insert(make_pair("Arg10", getParam_().getValue("labels:Arg10")));
    label_massshift_.insert(make_pair("Lys4", getParam_().getValue("labels:Lys4")));
    label_massshift_.insert(make_pair("Lys6", getParam_().getValue("labels:Lys6")));
    label_massshift_.insert(make_pair("Lys8", getParam_().getValue("labels:Lys8")));
    label_massshift_.insert(make_pair("Dimethyl0", getParam_().getValue("labels:Dimethyl0")));
    label_massshift_.insert(make_pair("Dimethyl4", getParam_().getValue("labels:Dimethyl4")));
    label_massshift_.insert(make_pair("Dimethyl6", getParam_().getValue("labels:Dimethyl6")));
    label_massshift_.insert(make_pair("Dimethyl8", getParam_().getValue("labels:Dimethyl8")));
    label_massshift_.insert(make_pair("ICPL0", getParam_().getValue("labels:ICPL0")));
    label_massshift_.insert(make_pair("ICPL4", getParam_().getValue("labels:ICPL4")));
    label_massshift_.insert(make_pair("ICPL6", getParam_().getValue("labels:ICPL6")));
    label_massshift_.insert(make_pair("ICPL10", getParam_().getValue("labels:ICPL10")));
  }

  /**
   * @brief generate list of mass patterns
   *
   * @return list of mass patterns
   */
  std::vector<MassPattern> generateMassPatterns_()
  {
    // SILAC, Dimethyl, ICPL or no labelling ??

    bool labelling_SILAC = ((labels_.find("Arg") != std::string::npos) || (labels_.find("Lys") != std::string::npos));
    bool labelling_Dimethyl = (labels_.find("Dimethyl") != std::string::npos);
    bool labelling_ICPL = (labels_.find("ICPL") != std::string::npos);
    bool labelling_none = labels_.empty() || (labels_ == "[]") || (labels_ == "()") || (labels_ == "{}");

    bool SILAC = (labelling_SILAC && !labelling_Dimethyl && !labelling_ICPL && !labelling_none);
    bool Dimethyl = (!labelling_SILAC && labelling_Dimethyl && !labelling_ICPL && !labelling_none);
    bool ICPL = (!labelling_SILAC && !labelling_Dimethyl && labelling_ICPL && !labelling_none);
    bool none = (!labelling_SILAC && !labelling_Dimethyl && !labelling_ICPL && labelling_none);

    if (!(SILAC || Dimethyl || ICPL || none))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Unknown labelling. Neither SILAC, Dimethyl nor ICPL.");
    }

    // debug output labels
    cout << "\n";
    for (unsigned i = 0; i < samples_labels_.size(); ++i)
    {
      cout << "sample " << (i + 1) << ":   ";
      for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
      {
        cout << samples_labels_[i][j] << " ";
      }
      cout << "\n";
    }

    // check if the labels are included in advanced section "labels"
    String all_labels = "Arg6 Arg10 Lys4 Lys6 Lys8 Dimethyl0 Dimethyl4 Dimethyl6 Dimethyl8 ICPL0 ICPL4 ICPL6 ICPL10";
    for (unsigned i = 0; i < samples_labels_.size(); i++)
    {
      for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
      {
        if (all_labels.find(samples_labels_[i][j]) == std::string::npos)
        {
          std::stringstream stream;
          stream << "The label " << samples_labels_[i][j] << " is unknown.";
          throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream.str());
        }
      }
    }

    // generate mass shift list
    std::vector<MassPattern> list;
    if (SILAC)
    {
      // SILAC
      // We assume the first sample to be unlabelled. Even if the "[]" for the first sample in the label string has not been specified.

      for (unsigned ArgPerPeptide = 0; ArgPerPeptide <= missed_cleavages_ + 1; ArgPerPeptide++)
      {
        for (unsigned LysPerPeptide = 0; LysPerPeptide <= missed_cleavages_ + 1; LysPerPeptide++)
        {
          if (ArgPerPeptide + LysPerPeptide <= missed_cleavages_ + 1)
          {
            MassPattern temp;
            temp.push_back(0);
            for (unsigned i = 0; i < samples_labels_.size(); i++)
            {
              double mass_shift = 0;
              // Considering the case of an amino acid (e.g. LysPerPeptide != 0) for which no label is present (e.g. Lys4There && Lys6There && Lys8There == false) makes no sense. Therefore each amino acid will have to give its "Go Ahead" before the shift is calculated.
              bool goAhead_Lys = false;
              bool goAhead_Arg = false;

              for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
              {
                bool Arg6There = (samples_labels_[i][j].find("Arg6") != std::string::npos); // Is Arg6 in the SILAC label?
                bool Arg10There = (samples_labels_[i][j].find("Arg10") != std::string::npos);
                bool Lys4There = (samples_labels_[i][j].find("Lys4") != std::string::npos);
                bool Lys6There = (samples_labels_[i][j].find("Lys6") != std::string::npos);
                bool Lys8There = (samples_labels_[i][j].find("Lys8") != std::string::npos);

                mass_shift = mass_shift + ArgPerPeptide * (Arg6There * label_massshift_["Arg6"] + Arg10There * label_massshift_["Arg10"]) + LysPerPeptide * (Lys4There * label_massshift_["Lys4"] + Lys6There * label_massshift_["Lys6"] + Lys8There * label_massshift_["Lys8"]);

                goAhead_Arg = goAhead_Arg || !(ArgPerPeptide != 0 && !Arg6There && !Arg10There);
                goAhead_Lys = goAhead_Lys || !(LysPerPeptide != 0 && !Lys4There && !Lys6There && !Lys8There);
              }

              if (goAhead_Arg && goAhead_Lys && (mass_shift != 0))
              {
                temp.push_back(mass_shift);
              }
            }

            if (temp.size() > 1)
            {
              list.push_back(temp);
            }
          }
        }
      }

    }
    else if (Dimethyl || ICPL)
    {
      // Dimethyl or ICPL
      // We assume each sample to be labelled only once.

      for (unsigned mc = 0; mc <= missed_cleavages_; ++mc)
      {
        MassPattern temp;
        for (unsigned i = 0; i < samples_labels_.size(); i++)
        {
          temp.push_back((mc + 1) * (label_massshift_[samples_labels_[i][0]] - label_massshift_[samples_labels_[0][0]]));
        }
        list.push_back(temp);
      }

    }
    else
    {
      // none (singlet detection)
      MassPattern temp;
      temp.push_back(0);
      list.push_back(temp);
    }

    // sort mass patterns
    // (from small mass shifts to larger ones, i.e. few miscleavages = simple explanation first)
    std::sort(list.begin(), list.end());

    // generate additional mass shifts due to knock-outs
    if (knock_out_ && list[0].size() == 1)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Knock-outs for singlet detection not relevant.");
    }
    else if (knock_out_ && list[0].size() <= 4)
    {
      generateKnockoutMassShifts(list);
    }
    else if (knock_out_ && list[0].size() > 4)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Knock-outs for multiplex experiments with more than 4 samples not supported.");
    }

    // debug output mass shifts
    cout << "\n";
    for (unsigned i = 0; i < list.size(); ++i)
    {
      std::cout << "mass shift " << (i + 1) << ":    ";
      MassPattern temp = list[i];
      for (unsigned j = 0; j < temp.size(); ++j)
      {
        std::cout << temp[j] << "  ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";

    return list;
  }

  /**
   * @brief generate all mass shifts that can occur due to the absence of one or multiple peptides
   * (e.g. for a triplet experiment generate the doublets and singlets that might be present)
   *
   * @param list of mass shifts to be extended
   */
  void generateKnockoutMassShifts(std::vector<MassPattern>& list)
  {
    unsigned n = list[0].size(); // n=2 for doublets, n=3 for triplets, n=4 for quadruplets
    if (knock_out_ && n == 4)
    {
      unsigned m = list.size();
      for (unsigned i = 0; i < m; ++i)
      {
        MassPattern triplet1(1, 0);
        triplet1.push_back(list[i][2] - list[i][1]);
        triplet1.push_back(list[i][3] - list[i][1]);
        list.push_back(triplet1);

        MassPattern triplet2(1, 0);
        triplet2.push_back(list[i][2] - list[i][0]);
        triplet2.push_back(list[i][3] - list[i][0]);
        list.push_back(triplet2);

        MassPattern triplet3(1, 0);
        triplet3.push_back(list[i][1] - list[i][0]);
        triplet3.push_back(list[i][2] - list[i][0]);
        list.push_back(triplet3);


        MassPattern doublet1(1, 0);
        doublet1.push_back(list[i][1]);
        list.push_back(doublet1);

        MassPattern doublet2(1, 0);
        doublet2.push_back(list[i][2]);
        list.push_back(doublet2);

        MassPattern doublet3(1, 0);
        doublet3.push_back(list[i][3]);
        list.push_back(doublet3);

        MassPattern doublet4(1, 0);
        doublet4.push_back(list[i][2] - list[i][1]);
        list.push_back(doublet4);

        MassPattern doublet5(1, 0);
        doublet5.push_back(list[i][3] - list[i][1]);
        list.push_back(doublet5);

        MassPattern doublet6(1, 0);
        doublet6.push_back(list[i][3] - list[i][2]);
        list.push_back(doublet6);
      }

      MassPattern singlet(1, 0);
      list.push_back(singlet);
    }
    else if (knock_out_ && n == 3)
    {
      unsigned m = list.size();
      for (unsigned i = 0; i < m; ++i)
      {
        MassPattern doublet1(1, 0);
        doublet1.push_back(list[i][1]);
        list.push_back(doublet1);

        MassPattern doublet2(1, 0);
        doublet2.push_back(list[i][2] - list[i][1]);
        list.push_back(doublet2);

        MassPattern doublet3(1, 0);
        doublet3.push_back(list[i][2]);
        list.push_back(doublet3);
      }

      MassPattern singlet(1, 0);
      list.push_back(singlet);
    }
    else if (knock_out_ && n == 2)
    {
      MassPattern singlet(1, 0);
      list.push_back(singlet);
    }
  }

  /**
   * @brief split labels string
   *
   * @return list of samples containing lists of corresponding labels
   */
  std::vector<std::vector<String> > splitLabelString_()
  {
    std::vector<std::vector<String> > samples_labels;
    std::vector<String> temp_samples;
    boost::split(temp_samples, labels_, boost::is_any_of("[](){}")); // any bracket allowed to separate samples
    for (unsigned i = 0; i < temp_samples.size(); ++i)
    {
      if (!temp_samples[i].empty())
      {
        vector<String> temp_labels;
        boost::split(temp_labels, temp_samples[i], boost::is_any_of(",;: ")); // various separators allowed to separate labels
        samples_labels.push_back(temp_labels);
      }
    }

    return samples_labels;
  }

  /**
   * @brief comparator of peak patterns
   *
   * @param pattern1    first peak pattern
   * @param pattern2    second peak pattern
   *
   * @return true if pattern1 should be searched before pattern2
   */
  static bool less_pattern(const MultiplexPeakPattern& pattern1, const MultiplexPeakPattern& pattern2)
  {
    if (pattern1.getMassShiftCount() == pattern2.getMassShiftCount())
    {
      if (pattern1.getCharge() == pattern2.getCharge())
      {
        // The first mass shift is by definition always zero.
        if ((pattern1.getMassShiftCount() > 1) && (pattern2.getMassShiftCount() > 1))
        {
          // 4Da before 8Da etc. (larger miss cleavages last)
          return pattern1.getMassShiftAt(1) < pattern2.getMassShiftAt(1);
        }
        else
        {
          // Should never happen.
          return true;
        }
      }
      else
      {
        // 5+ before 4+ before 3+ etc.
        return pattern1.getCharge() > pattern2.getCharge();
      }
    }
    else
    {
      // triplets before doublets before singlets
      return pattern1.getMassShiftCount() > pattern2.getMassShiftCount();
    }
  }

  /**
   * @brief generate list of mass shifts
   *
   * @param charge_min    minimum charge
   * @param charge_max    maximum charge
   * @param peaks_per_peptide_max    maximum number of isotopes in peptide
   * @param mass_pattern_list    mass shifts due to labelling
   *
   * @return list of mass shifts
   */
  std::vector<MultiplexPeakPattern> generatePeakPatterns_(int charge_min, int charge_max, int peaks_per_peptide_max, std::vector<MassPattern> mass_pattern_list)
  {
    std::vector<MultiplexPeakPattern> list;

    // iterate over all charge states
    for (int c = charge_max; c >= charge_min; --c)
    {
      // iterate over all mass shifts
      for (unsigned i = 0; i < mass_pattern_list.size(); ++i)
      {
        MultiplexPeakPattern pattern(c, peaks_per_peptide_max, mass_pattern_list[i], i);
        list.push_back(pattern);
      }
    }
    
    sort(list.begin(),list.end(),less_pattern);
    
    // debug output
    /*for (size_t i = 0; i < list.size(); ++i)
    {
      std::cout << list[i].getCharge() << "+  ";
      for (size_t j = 0; j < list[i].getMassShiftCount(); ++j)
      {
        std::cout << list[i].getMassShiftAt(j) << "  ";
      }
      std::cout << "\n";
    }*/

    return list;
  }

  /**
   * @brief calculate peptide intensities
   *
   * @param all_intensities    vectors of profile intensities for each of the peptides (first index: peptide 0=L, 1=M, 2=H etc, second index: raw data point)
   * @return vector with intensities for each peptide
   */
  std::vector<double> getPeptideIntensities(std::vector<std::vector<double> >& all_intensities)
  {
    OPENMS_PRECONDITION(!all_intensities.empty(), "The entire profile intensity vector should not be empty.");
    bool empty_intensities = false;
    for (unsigned i = 0; i < all_intensities.size(); ++i)
    {
      empty_intensities = empty_intensities || all_intensities[i].empty();
    }
    OPENMS_PRECONDITION(!empty_intensities, "None of the individual profile intensity vectors should be empty.");
    unsigned count = all_intensities[0].size();
    for (unsigned i = 0; i < all_intensities.size(); ++i)
    {
      if (all_intensities[i].size() != count)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "The profile intensity vectors for each peptide are not of the same size.");
      }
    }

    // determine ratios through linear regression
    // of all (spline-interpolated) profile intensities that passed all multiplex filters, @see MultiplexFiltering
    std::vector<double> ratios; // L:L, M:L, H:L etc.
    std::vector<double> intensities; // L, M, H etc.
    // loop over peptides
    for (unsigned i = 0; i < all_intensities.size(); ++i)
    {
      // filter for non-NaN intensities
      std::vector<double> intensities1;
      std::vector<double> intensities2;
      double intensity = 0;
      for (unsigned j = 0; j < count; ++j)
      {
        if (!(boost::math::isnan(all_intensities[0][j])) && !(boost::math::isnan(all_intensities[i][j])))
        {
          intensities1.push_back(all_intensities[0][j]);
          intensities2.push_back(all_intensities[i][j]);
          intensity += all_intensities[i][j];
        }
      }

      LinearRegressionWithoutIntercept linreg;
      linreg.addData(intensities1, intensities2);
      ratios.push_back(linreg.getSlope());
      intensities.push_back(intensity);

    }

    // correct peptide intensities
    // The peptide ratios are calculated as linear regression of (spline-interpolated) profile intensities, @see linreg
    // The individual peptide intensities are the sum of the same profile intensities. But the quotient of these peptide intensities
    // is not necessarily the same as the independently calculated ratio from the linear regression. Since the peptide ratio
    // from linear regression is the more accurate one, we correct the two peptide intensities by projecting them onto the ratio.
    // In the end, both peptide ratio from linear regression and the quotient of the peptide intensities are identical.
    std::vector<double> corrected_intensities;
    if (all_intensities.size() == 2)
    {
      double intensity1 = (intensities[0] + ratios[1] * intensities[1]) / (1 + ratios[1] * ratios[1]);
      double intensity2 = ratios[1] * intensity1;
      corrected_intensities.push_back(intensity1);
      corrected_intensities.push_back(intensity2);
    }
    else if (all_intensities.size() > 2)
    {
      // Now with n instead of two peptide intensities, one needs to project the peptide intensities onto the hyperplane defined
      // by the set of all peptide ratios (TODO). Instead, it is simpler to keep the lightest peptide intensity fixed, and correct
      // only the remaining ones. The correct peptide ratio (from linear regression) is reported on both cases.
      corrected_intensities.push_back(intensities[0]);
      for (unsigned i = 1; i < all_intensities.size(); ++i)
      {
        corrected_intensities.push_back(ratios[i] * intensities[0]);
      }
    }
    else
    {
      // For simple feature detection (singlets) the intensities remain unchanged.
      corrected_intensities.push_back(intensities[0]);
    }

    return corrected_intensities;
  }

  /**
   * @brief generates consensus and feature maps containing all peptide multiplets
   *
   * @param centroided    type of spectral input data (profile or centroided)
   * @param patterns    patterns of isotopic peaks we have been searching for
   * @param filter_results    filter results for each of the patterns
   * @param cluster_results    clusters of filter results
   * @param consensus_map    consensus map with peptide multiplets (to be filled)
   * @param feature_map    feature map with peptides (to be filled)
   */
  void generateMaps_(bool centroided, std::vector<MultiplexPeakPattern> patterns, std::vector<MultiplexFilterResult> filter_results, std::vector<std::map<int, GridBasedCluster> > cluster_results, ConsensusMap& consensus_map, FeatureMap& feature_map)
  {
    // loop over peak patterns
    for (unsigned pattern = 0; pattern < patterns.size(); ++pattern)
    {
      // loop over clusters
      for (std::map<int, GridBasedCluster>::const_iterator cluster_it = cluster_results[pattern].begin(); cluster_it != cluster_results[pattern].end(); ++cluster_it)
      {
        ConsensusFeature consensus;

        // The position (m/z, RT) of the peptide features is the centre-of-mass of the mass trace of the lightest isotope.
        // The centre-of-mass is the intensity-weighted average of the peak positions.
        unsigned number_of_peptides = patterns[pattern].getMassShiftCount();
        std::vector<double> sum_intensity_mz(number_of_peptides, 0);
        std::vector<double> sum_intensity_rt(number_of_peptides, 0);
        std::vector<double> sum_intensity(number_of_peptides, 0);
        // intensities for ratio determination.
        // For centroided input data, these intensities are peak intensities.
        // For profile input data, these intensities are the spline-interpolated profile intensities.
        // First index is the peptide, second is just the profile intensities collected
        std::vector<std::vector<double> > all_intensities(patterns[pattern].getMassShiftCount(), std::vector<double>());
        // bounding boxes of mass traces for each peptide multiplet
        // First index is the peptide, second is the mass trace within the peptide.
        std::map<std::pair<unsigned, unsigned>, DBoundingBox<2> > mass_traces;

        GridBasedCluster cluster = cluster_it->second;
        std::vector<int> points = cluster.getPoints();

        // loop over points in cluster
        for (std::vector<int>::const_iterator point_it = points.begin(); point_it != points.end(); ++point_it)
        {
          int index = (*point_it);

          MultiplexFilterResultPeak result_peak = filter_results[pattern].getFilterResultPeak(index);
          double rt = result_peak.getRT();

          for (unsigned peptide = 0; peptide < patterns[pattern].getMassShiftCount(); ++peptide)
          {
            sum_intensity_mz[peptide] += (result_peak.getMZ() + result_peak.getMZShifts()[(isotopes_per_peptide_max_ + 1) * peptide + 1]) * result_peak.getIntensities()[(isotopes_per_peptide_max_ + 1) * peptide + 1];
            sum_intensity_rt[peptide] += result_peak.getRT() * result_peak.getIntensities()[(isotopes_per_peptide_max_ + 1) * peptide + 1];
            sum_intensity[peptide]  += result_peak.getIntensities()[(isotopes_per_peptide_max_ + 1) * peptide + 1];
          }

          if (centroided)
          {
            // loop over isotopic peaks in peptide
            for (unsigned peak = 0; peak < isotopes_per_peptide_max_; ++peak)
            {
              // loop over peptides
              for (unsigned peptide = 0; peptide < patterns[pattern].getMassShiftCount(); ++peptide)
              {
                unsigned index = (isotopes_per_peptide_max_ + 1) * peptide + peak + 1; // +1 due to zeroth peaks
                all_intensities[peptide].push_back(result_peak.getIntensities()[index]); // Note that the intensity can be NaN. To be checked later.

                double mz_shift = result_peak.getMZShifts()[index];
                if (!(boost::math::isnan(mz_shift)))
                {
                  std::pair<unsigned, unsigned> peptide_peak(peptide, peak);
                  mass_traces[peptide_peak].enlarge(rt, result_peak.getMZ() + mz_shift);
                }
              }
            }
          }
          else
          {
            // iterate over profile data
            // (We use the (spline-interpolated) profile intensities for a very accurate ratio determination.)
            for (int i = 0; i < result_peak.size(); ++i)
            {
              MultiplexFilterResultRaw result_raw = result_peak.getFilterResultRaw(i);

              // loop over isotopic peaks in peptide
              for (unsigned peak = 0; peak < isotopes_per_peptide_max_; ++peak)
              {
                // loop over peptides
                for (unsigned peptide = 0; peptide < patterns[pattern].getMassShiftCount(); ++peptide)
                {
                  unsigned index = (isotopes_per_peptide_max_ + 1) * peptide + peak + 1; // +1 due to zeroth peaks
                  all_intensities[peptide].push_back(result_raw.getIntensities()[index]); // Note that the intensity can be NaN. To be checked later.

                  double mz_shift = result_raw.getMZShifts()[index];
                  if (!(boost::math::isnan(mz_shift)))
                  {
                    std::pair<unsigned, unsigned> peptide_peak(peptide, peak);
                    mass_traces[peptide_peak].enlarge(rt, result_raw.getMZ() + mz_shift);
                  }
                }
              }
            }
          }

        }

        // calculate intensities for each of the peptides from profile data
        std::vector<double> peptide_intensities = getPeptideIntensities(all_intensities);

        // average peptide intensity (= consensus intensity)
        double average_peptide_intensity = 0;
        for (unsigned i = 0; i < peptide_intensities.size(); ++i)
        {
          average_peptide_intensity += peptide_intensities[i];
        }
        average_peptide_intensity /= peptide_intensities.size();

        // fill map with consensuses and its features
        consensus.setMZ(sum_intensity_mz[0] / sum_intensity[0]);
        consensus.setRT(sum_intensity_rt[0] / sum_intensity[0]);
        consensus.setIntensity(average_peptide_intensity);
        consensus.setCharge(patterns[pattern].getCharge());
        consensus.setQuality(1 - 1 / points.size()); // rough quality score in [0,1]

        for (unsigned peptide = 0; peptide < patterns[pattern].getMassShiftCount(); ++peptide)
        {
          FeatureHandle feature_handle;
          feature_handle.setMZ(sum_intensity_mz[peptide] / sum_intensity[peptide]);
          feature_handle.setRT(sum_intensity_rt[peptide] / sum_intensity[peptide]);
          feature_handle.setIntensity(peptide_intensities[peptide]);
          feature_handle.setCharge(patterns[pattern].getCharge());
          feature_handle.setMapIndex(peptide);
          //feature_handle.setUniqueId(&UniqueIdInterface::setUniqueId);    // TODO: Do we need to set unique ID?
          consensus_map.getFileDescriptions()[peptide].size++;
          consensus.insert(feature_handle);

          Feature feature;
          feature.setMZ(sum_intensity_mz[peptide] / sum_intensity[peptide]);
          feature.setRT(sum_intensity_rt[peptide] / sum_intensity[peptide]);
          feature.setIntensity(peptide_intensities[peptide]);
          feature.setCharge(patterns[pattern].getCharge());
          feature.setOverallQuality(1 - 1 / points.size());
          for (unsigned peak = 0; peak < isotopes_per_peptide_max_; ++peak)
          {
            std::pair<unsigned, unsigned> peptide_peak(peptide, peak);
            if (mass_traces.count(peptide_peak) > 0)
            {
              ConvexHull2D hull;
              hull.addPoint(DPosition<2>(mass_traces[peptide_peak].minX(), mass_traces[peptide_peak].minY()));
              hull.addPoint(DPosition<2>(mass_traces[peptide_peak].minX(), mass_traces[peptide_peak].maxY()));
              hull.addPoint(DPosition<2>(mass_traces[peptide_peak].maxX(), mass_traces[peptide_peak].minY()));
              hull.addPoint(DPosition<2>(mass_traces[peptide_peak].maxX(), mass_traces[peptide_peak].maxY()));
              feature.getConvexHulls().push_back(hull);
            }
          }

          feature_map.push_back(feature);
        }

        consensus_map.push_back(consensus);
      }

    }

  }

  /**
   * @brief generates the data structure for mzQuantML output
   *
   * @param consensus_map    consensus map with complete quantitative information
   * @param quantifications    MSQuantifications data structure for writing mzQuantML (mzq)
   */
  void generateMSQuantifications(MSExperiment<Peak1D>& exp, ConsensusMap& consensus_map, MSQuantifications& quantifications)
  {
    // generate the labels
    // (for each sample a list of (label string, mass shift) pairs)
    // for example triple-SILAC: [(none,0)][(Lys4,4.0251),(Arg6,6.0201)][Lys8,8.0141)(Arg10,10.0082)]
    std::vector<std::vector<std::pair<String, double> > > labels;
    for (unsigned sample = 0; sample < samples_labels_.size(); ++sample)
    {
      // The labels are required to be ordered in mass shift.
      std::map<double, String> single_label_map;
      std::vector<std::pair<String, double> > single_label;
      for (unsigned label = 0; label < samples_labels_[sample].size(); ++label)
      {
        String label_string = samples_labels_[sample][label];
        double shift;
        if (label_string == "")
        {
          label_string = "none";
          shift = 0;
        }
        else
        {
          shift = label_massshift_[label_string];
        }

        single_label_map[shift] = label_string;
      }
      for (std::map<double, String>::const_iterator it = single_label_map.begin(); it != single_label_map.end(); ++it)
      {
        std::pair<String, double> label_shift(it->second, it->first);
        single_label.push_back(label_shift);
      }
      labels.push_back(single_label);
    }

    quantifications.registerExperiment(exp, labels);
    quantifications.assignUIDs();

    MSQuantifications::QUANT_TYPES quant_type = MSQuantifications::MS1LABEL;
    quantifications.setAnalysisSummaryQuantType(quant_type);

    // add results from  analysis
    LOG_DEBUG << "Generating output mzQuantML file..." << endl;
    ConsensusMap numap(consensus_map);
    //calculate ratios
    for (ConsensusMap::iterator cit = numap.begin(); cit != numap.end(); ++cit)
    {
      // make ratio templates
      std::vector<ConsensusFeature::Ratio> rts;
      for (std::vector<MSQuantifications::Assay>::const_iterator ait = quantifications.getAssays().begin() + 1; ait != quantifications.getAssays().end(); ++ait)
      {
        ConsensusFeature::Ratio r;
        r.numerator_ref_ = String(quantifications.getAssays().begin()->uid_);
        r.denominator_ref_ = String(ait->uid_);
        r.description_.push_back("Simple ratio calc");
        r.description_.push_back("light to medium/.../heavy");
        rts.push_back(r);
      }

      const ConsensusFeature::HandleSetType& feature_handles = cit->getFeatures();
      if (feature_handles.size() > 1)
      {
        std::set<FeatureHandle, FeatureHandle::IndexLess>::const_iterator fit = feature_handles.begin(); // this is unlabeled
        ++fit;
        for (; fit != feature_handles.end(); ++fit)
        {
          Size ri = std::distance(feature_handles.begin(), fit);
          rts[ri - 1].ratio_value_ =  feature_handles.begin()->getIntensity() / fit->getIntensity(); // a proper algo should never have 0-intensities so no 0devison ...
        }
      }

      cit->setRatios(rts);
    }
    quantifications.addConsensusMap(numap); //add FeatureFinderMultiplex result

  }

  /**
   * @brief Write consensus map to consensusXML file.
   *
   * @param filename    name of consensusXML file
   * @param map    consensus map for output
   */
  void writeConsensusMap_(const String& filename, ConsensusMap& map) const
  {
    map.sortByPosition();
    map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    map.setExperimentType("multiplex");

    // annotate maps
    for (unsigned i = 0; i < samples_labels_.size(); ++i)
    {
      ConsensusMap::FileDescription& desc = map.getFileDescriptions()[i];
      desc.filename = filename;

      if (knock_out_)
      {
        // With knock-outs present, the correct labels can only be determined during ID mapping.
        desc.label = "";
      }
      else
      {
        String label_string;
        for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
        {
          label_string.append(samples_labels_[i][j]);
        }
        desc.label = label_string;
      }
    }

    ConsensusXMLFile file;
    file.store(filename, map);
  }

  /**
   * @brief Write feature map to featureXML file.
   *
   * @param filename    name of featureXML file
   * @param map    feature map for output
   */
  void writeFeatureMap_(const String& filename, FeatureMap& map) const
  {
    map.sortByPosition();
    map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    FeatureXMLFile file;
    file.store(filename, map);
  }

  /**
   * @brief Write MS quantification map to mzq file.
   *
   * @param filename    name of mzq file
   * @param map    MS quantification map for output
   */
  void writeMSQuantifications(const String& filename, MSQuantifications& msq) const
  {
    MzQuantMLFile file;
    file.store(filename, msq);
  }

  /**
  * @brief simple linear regression through the origin
  *
  * TODO: combine with OpenMS/MATH/STATISTICS/LinearRegression
  */
  class LinearRegressionWithoutIntercept
  {
public:
    /**
     * @brief constructor
     */
    LinearRegressionWithoutIntercept() :
      sum_xx_(0), sum_xy_(0), n_(0)
    {
    }

    /**
     * @brief adds an observation (x,y) to the regression data set.
     *
     * @param x    independent variable value
     * @param y    dependent variable value
     */
    void addData(double x, double y)
    {
      sum_xx_ += x * x;
      sum_xy_ += x * y;

      ++n_;
    }

    /**
     * @brief adds observations (x,y) to the regression data set.
     *
     * @param x    vector of independent variable values
     * @param y    vector of dependent variable values
     */
    void addData(std::vector<double>& x, std::vector<double>& y)
    {
      for (unsigned i = 0; i < x.size(); ++i)
      {
        addData(x[i], y[i]);
      }
    }

    /**
     * @brief returns the slope of the estimated regression line.
     */
    double getSlope()
    {
      if (n_ < 2)
      {
        return std::numeric_limits<double>::quiet_NaN(); // not enough data
      }
      return sum_xy_ / sum_xx_;
    }

private:
    /**
     * @brief total variation in x
     */
    double sum_xx_;

    /**
     * @brief sum of products
     */
    double sum_xy_;

    /**
     * @brief number of observations
     */
    int n_;

  };


  ExitCodes main_(int, const char**)
  {

    /**
     * handle parameters
     */
    getParameters_in_out_();
    getParameters_labels_();
    getParameters_algorithm_();

    /**
     * load input
     */
    MzMLFile file;
    MSExperiment<Peak1D> exp;

    // only read MS1 spectra
    std::vector<int> levels;
    levels.push_back(1);
    file.getOptions().setMSLevels(levels);

    LOG_DEBUG << "Loading input..." << endl;
    file.setLogType(log_type_);
    file.load(in_, exp);

    // update m/z and RT ranges
    exp.updateRanges();

    // sort according to RT and MZ
    exp.sortSpectra();

    // determine type of spectral data (profile or centroided)
    bool centroided = (PeakTypeEstimator().estimateType(exp[0].begin(), exp[0].end()) == SpectrumSettings::PEAKS);

    /**
     * pick peaks
     */
    MSExperiment<Peak1D> exp_picked;
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s; // peak boundaries for spectra
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c; // peak boundaries for chromatograms
    if (!centroided)
    {
      PeakPickerHiRes picker;
      Param param = picker.getParameters();
      picker.setLogType(log_type_);
      param.setValue("ms_levels", ListUtils::create<Int>("1"));
      param.setValue("signal_to_noise", 0.0); // signal-to-noise estimation switched off
      picker.setParameters(param);

      picker.pickExperiment(exp, exp_picked, boundaries_exp_s, boundaries_exp_c);
    }

    /**
     * filter for peak patterns
     */
    bool missing_peaks_ = false;
    std::vector<MassPattern> masses = generateMassPatterns_();
    std::vector<MultiplexPeakPattern> patterns = generatePeakPatterns_(charge_min_, charge_max_, isotopes_per_peptide_max_, masses);

    std::vector<MultiplexFilterResult> filter_results;
    if (centroided)
    {
      // centroided data
      MultiplexFilteringCentroided filtering(exp, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, missing_peaks_, intensity_cutoff_, mz_tolerance_, mz_unit_, peptide_similarity_, averagine_similarity_, averagine_similarity_scaling_, out_debug_);
      filtering.setLogType(log_type_);
      filter_results = filtering.filter();
    }
    else
    {
      // profile data
      MultiplexFilteringProfile filtering(exp, exp_picked, boundaries_exp_s, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, missing_peaks_, intensity_cutoff_, mz_tolerance_, mz_unit_, peptide_similarity_, averagine_similarity_, averagine_similarity_scaling_, out_debug_);
      filtering.setLogType(log_type_);
      filter_results = filtering.filter();
    }

    /**
     * cluster filter results
     */
    std::vector<std::map<int, GridBasedCluster> > cluster_results;
    if (centroided)
    {
      // centroided data
      MultiplexClustering clustering(exp, mz_tolerance_, mz_unit_, rt_typical_, rt_min_, out_debug_);
      clustering.setLogType(log_type_);
      cluster_results = clustering.cluster(filter_results);
    }
    else
    {
      // profile data
      MultiplexClustering clustering(exp, exp_picked, boundaries_exp_s, rt_typical_, rt_min_, out_debug_);
      clustering.setLogType(log_type_);
      cluster_results = clustering.cluster(filter_results);
    }

    /**
     * write to output
     */
    ConsensusMap consensus_map;
    FeatureMap feature_map;
    generateMaps_(centroided, patterns, filter_results, cluster_results, consensus_map, feature_map);
    if (out_ != "")
    {
      writeConsensusMap_(out_, consensus_map);
    }
    if (out_features_ != "")
    {
      writeFeatureMap_(out_features_, feature_map);
    }

    if (out_mzq_ != "")
    {
      MSQuantifications quantifications;
      generateMSQuantifications(exp, consensus_map, quantifications);
      writeMSQuantifications(out_mzq_, quantifications);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderMultiplex tool;
  return tool.main(argc, argv);
}

//@endcond
