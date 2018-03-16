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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Timo Sachsenberg, Samuel Wein, Mathias Walzer $
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/COMPARISON/CLUSTERING/GridBasedCluster.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

//Contrib includes
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
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

  @brief Detects peptide pairs in LC-MS data and determines their relative abundance.

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

  FeatureFinderMultiplex is a tool for the fully automated analysis of quantitative proteomics data. It detects pairs of isotopic envelopes with fixed m/z separation. It requires no prior sequence identification of the peptides. In what follows we outline the algorithm.

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
  String spectrum_type_;
  String averagine_type_;

  // section "labels"
  map<String, double> label_mass_shift_;

public:
  TOPPFeatureFinderMultiplex() :
    TOPPBase("FeatureFinderMultiplex", "Determination of peak ratios in LC-MS data", true),
    charge_min_(1), charge_max_(1), missed_cleavages_(0), isotopes_per_peptide_min_(1), isotopes_per_peptide_max_(1), rt_typical_(0.0), rt_min_(0.0),
    mz_tolerance_(0.0), mz_unit_(true), intensity_cutoff_(0.0), peptide_similarity_(0.0), averagine_similarity_(0.0), averagine_similarity_scaling_(0.0), knock_out_(false)
  {
  }

  typedef std::vector<double> MassPattern; // list of mass shifts

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "LC-MS dataset in centroid or profile mode");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Set of all identified peptide groups (i.e. peptide pairs or triplets or singlets or ..). The m/z-RT positions correspond to the lightest peptide in each group.", false);
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out_features", "<file>", "", "Optional output file containing the individual peptide features in \'out\'.", false, true);
    setValidFormats_("out_features", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out_mzq", "<file>", "", "Optional output file of MzQuantML.", false, true);
    setValidFormats_("out_mzq", ListUtils::create<String>("mzq"));

    registerSubsection_("algorithm", "Parameters for the algorithm.");
    registerSubsection_("labels", "Isotopic labels that can be specified in section \'algorithm:labels\'.");

  }

  // create parameters for sections (set default values and restrictions)
  Param getSubsectionDefaults_(const String& section) const override
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
      defaults.setMinFloat("peptide_similarity", -1.0);
      defaults.setMaxFloat("peptide_similarity", 1.0);
      defaults.setValue("averagine_similarity", 0.4, "The isotopic pattern of a peptide should resemble the averagine model at this m/z position. This parameter is a lower bound on similarity between measured isotopic pattern and the averagine model.");
      defaults.setMinFloat("averagine_similarity", -1.0);
      defaults.setMaxFloat("averagine_similarity", 1.0);
      defaults.setValue("averagine_similarity_scaling", 0.75, "Let x denote this scaling factor, and p the averagine similarity parameter. For the detection of single peptides, the averagine parameter p is replaced by p' = p + x(1-p), i.e. x = 0 -> p' = p and x = 1 -> p' = 1. (For knock_out = true, peptide doublets and singlets are detected simulataneously. For singlets, the peptide similarity filter is irreleavant. In order to compensate for this 'missing filter', the averagine parameter p is replaced by the more restrictive p' when searching for singlets.)", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("averagine_similarity_scaling", 0.0);
      defaults.setMaxFloat("averagine_similarity_scaling", 1.0);
      defaults.setValue("missed_cleavages", 0, "Maximum number of missed cleavages due to incomplete digestion. (Only relevant if enzymatic cutting site coincides with labelling site. For example, Arg/Lys in the case of trypsin digestion and SILAC labelling.)");
      defaults.setMinInt("missed_cleavages", 0);
      defaults.setValue("knock_out", "false", "Is it likely that knock-outs are present? (Supported for doublex, triplex and quadruplex experiments only.)", ListUtils::create<String>("advanced"));
      defaults.setValidStrings("knock_out", ListUtils::create<String>("true,false"));
      defaults.setValue("spectrum_type", "automatic", "Type of MS1 spectra in input mzML file. 'automatic' determines the spectrum type directly from the input mzML file.", ListUtils::create<String>("advanced"));
      defaults.setValidStrings("spectrum_type", ListUtils::create<String>("profile,centroid,automatic"));
      defaults.setValue("averagine_type","peptide","The type of averagine to use, currently RNA, DNA or peptide", ListUtils::create<String>("advanced"));
      defaults.setValidStrings("averagine_type", ListUtils::create<String>("peptide,RNA,DNA"));
    }

    if (section == "labels")
    {
      MultiplexDeltaMassesGenerator generator;
      Param p = generator.getParameters();
      
      for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
      {
        defaults.setValue(it->name, it->value, it->description, ListUtils::create<String>("advanced"));
        defaults.setMinFloat(it->name, 0.0);
      }
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
    spectrum_type_ = getParam_().getValue("algorithm:spectrum_type");
    averagine_type_ = getParam_().getValue("algorithm:averagine_type");
  }

  /**
   * @brief process parameters of 'labels' section
   */
  void getParameters_labels_()
  {
    Param p = getParam_();
    
    // create map of pairs (label as string, mass shift as double)
    for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
    {
      label_mass_shift_.insert(make_pair(it->name, it->value));
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
    
    String labels(labels_);
    boost::replace_all(labels, "[]", "no_label");
    boost::replace_all(labels, "()", "no_label");
    boost::replace_all(labels, "{}", "no_label");
    boost::split(temp_samples, labels, boost::is_any_of("[](){}")); // any bracket allowed to separate samples
    
    for (unsigned i = 0; i < temp_samples.size(); ++i)
    {
      if (!temp_samples[i].empty())
      {
        if (temp_samples[i]=="no_label")
        {
          vector<String> temp_labels;
          temp_labels.push_back("no_label");
          samples_labels.push_back(temp_labels);
        }
        else
        {
          vector<String> temp_labels;
          boost::split(temp_labels, temp_samples[i], boost::is_any_of(",;: ")); // various separators allowed to separate labels
          samples_labels.push_back(temp_labels);
        }
      }
    }
    
    if (samples_labels.empty())
    {
      vector<String> temp_labels;
      temp_labels.push_back("no_label");
      samples_labels.push_back(temp_labels);
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
  static bool less_pattern(const MultiplexIsotopicPeakPattern& pattern1, const MultiplexIsotopicPeakPattern& pattern2)
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
   * @brief generate list of m/z shifts
   *
   * @param charge_min    minimum charge
   * @param charge_max    maximum charge
   * @param peaks_per_peptide_max    maximum number of isotopes in peptide
   * @param mass_pattern_list    mass shifts due to labelling
   *
   * @return list of m/z shifts
   */
  std::vector<MultiplexIsotopicPeakPattern> generatePeakPatterns_(int charge_min, int charge_max, int peaks_per_peptide_max, std::vector<MultiplexDeltaMasses> mass_pattern_list)
  {
    std::vector<MultiplexIsotopicPeakPattern> list;

    // iterate over all charge states
    for (int c = charge_max; c >= charge_min; --c)
    {
      // iterate over all mass shifts
      for (unsigned i = 0; i < mass_pattern_list.size(); ++i)
      {
        MultiplexIsotopicPeakPattern pattern(c, peaks_per_peptide_max, mass_pattern_list[i], i);
        list.push_back(pattern);
      }
    }
    
    sort(list.begin(),list.end(),less_pattern);
    
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
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The profile intensity vectors for each peptide are not of the same size.");
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
  void generateMaps_(bool centroided, std::vector<MultiplexIsotopicPeakPattern> patterns, std::vector<MultiplexFilterResult> filter_results, std::vector<std::map<int, GridBasedCluster> > cluster_results, ConsensusMap& consensus_map, FeatureMap& feature_map)
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
          int idx = (*point_it);

          MultiplexFilterResultPeak result_peak = filter_results[pattern].getFilterResultPeak(idx);
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
   * @param exp    experimental data
   * @param consensus_map    consensus map with complete quantitative information
   * @param quantifications    MSQuantifications data structure for writing mzQuantML (mzq)
   */
  void generateMSQuantifications(PeakMap& exp, ConsensusMap& consensus_map, MSQuantifications& quantifications)
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
          shift = label_mass_shift_[label_string];
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
    map.setExperimentType("labeled_MS1");

    // annotate maps
    for (unsigned i = 0; i < samples_labels_.size(); ++i)
    {
      ConsensusMap::FileDescription& desc = map.getFileDescriptions()[i];
      desc.filename = filename;

      if (knock_out_)
      {
        // With knock-outs present, the correct labels can only be determined during ID mapping.
        // For now, we simply store a unique identifier.
        std::stringstream stream;
        stream << "label " << i;
        desc.label = stream.str();
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


  ExitCodes main_(int, const char**) override
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
    PeakMap exp;

    // only read MS1 spectra
    std::vector<int> levels;
    levels.push_back(1);
    file.getOptions().setMSLevels(levels);

    LOG_DEBUG << "Loading input..." << endl;
    file.setLogType(log_type_);
    file.load(in_, exp);

    if (exp.getSpectra().empty())
    {
      throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS1 spectra in input file.");
    }

    // update m/z and RT ranges
    exp.updateRanges();

    // sort according to RT and MZ
    exp.sortSpectra();

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = exp[0].getType();
    if (spectrum_type == SpectrumSettings::UNKNOWN)
    {
      spectrum_type = PeakTypeEstimator().estimateType(exp[0].begin(), exp[0].end());
    }

    bool centroided;
    if (spectrum_type_=="automatic")
    {
      centroided = spectrum_type == SpectrumSettings::CENTROID;
    }
    else if (spectrum_type_=="centroid")
    {
      centroided = true;
    }
    else  // "profile"
    {
      centroided = false;
    }

    /**
     * pick peaks
     */
    PeakMap exp_picked;
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
    MultiplexDeltaMassesGenerator generator = MultiplexDeltaMassesGenerator(labels_, missed_cleavages_, label_mass_shift_);
    if (knock_out_)
    {
      generator.generateKnockoutDeltaMasses();
    }
    generator.printSamplesLabelsList();
    generator.printDeltaMassesList();
    
    std::vector<MultiplexDeltaMasses> masses = generator.getDeltaMassesList();
    std::vector<MultiplexIsotopicPeakPattern> patterns = generatePeakPatterns_(charge_min_, charge_max_, isotopes_per_peptide_max_, masses);

    bool missing_peaks_ = false;
    std::vector<MultiplexFilterResult> filter_results;
    if (centroided)
    {
      // centroided data
      MultiplexFilteringCentroided filtering(exp, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, missing_peaks_, intensity_cutoff_, mz_tolerance_, mz_unit_, peptide_similarity_, averagine_similarity_, averagine_similarity_scaling_, averagine_type_);
      filtering.setLogType(log_type_);
      filter_results = filtering.filter();
    }
    else
    {
      // profile data
      MultiplexFilteringProfile filtering(exp, exp_picked, boundaries_exp_s, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, missing_peaks_, intensity_cutoff_, mz_tolerance_, mz_unit_, peptide_similarity_, averagine_similarity_, averagine_similarity_scaling_, averagine_type_);
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
      MultiplexClustering clustering(exp, mz_tolerance_, mz_unit_, rt_typical_, rt_min_);
      clustering.setLogType(log_type_);
      cluster_results = clustering.cluster(filter_results);
    }
    else
    {
      // profile data
      MultiplexClustering clustering(exp, exp_picked, boundaries_exp_s, rt_typical_, rt_min_);
      clustering.setLogType(log_type_);
      cluster_results = clustering.cluster(filter_results);
    }

    /**
     * write to output
     */
    ConsensusMap consensus_map;
    StringList ms_runs;
    exp.getPrimaryMSRunPath(ms_runs);
    consensus_map.setPrimaryMSRunPath(ms_runs);

    FeatureMap feature_map;
    feature_map.setPrimaryMSRunPath(ms_runs);

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
