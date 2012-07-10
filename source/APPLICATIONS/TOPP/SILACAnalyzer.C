// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

class TOPPSILACAnalyzer
: public TOPPBase
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

    // section "labels"
    map<String, DoubleReal> label_identifiers;
    vector<vector <String> > SILAClabels;     // list of SILAC labels, e.g. selected_labels="[Lys4,Arg6][Lys8,Arg10]" => SILAClabels[0][1]="Arg6"
    vector<vector <DoubleReal> > massShifts;      // list of mass shifts

    typedef SILACClustering Clustering;

    vector<vector<SILACPattern> > data;
    vector<Clustering *> cluster_data;

	MSQuantifications msq;

  public:
    TOPPSILACAnalyzer()
      : TOPPBase("SILACAnalyzer","Determination of peak ratios in LC-MS data",true), allow_missing_peaks(true)
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
  Param getSubsectionDefaults_(const String& section) const
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
      defaults.setValue("missed_cleavages", 0 , "Maximum number of missed cleavages.");
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

    // create iterators for all labels to get corresponding mass shift
    map<String,DoubleReal>::iterator arg6 = label_identifiers.find("Arg6");
    map<String,DoubleReal>::iterator arg10 = label_identifiers.find("Arg10");
    map<String,DoubleReal>::iterator lys4 = label_identifiers.find("Lys4");
    map<String,DoubleReal>::iterator lys6 = label_identifiers.find("Lys6");
    map<String,DoubleReal>::iterator lys8 = label_identifiers.find("Lys8");
    map<String,DoubleReal>::iterator methyl4 = label_identifiers.find("Methyl4");
    map<String,DoubleReal>::iterator methyl8 = label_identifiers.find("Methyl8");
    map<String,DoubleReal>::iterator methyl12 = label_identifiers.find("Methyl12");
    map<String,DoubleReal>::iterator methyl16 = label_identifiers.find("Methyl16");
    map<String,DoubleReal>::iterator methyl24 = label_identifiers.find("Methyl24");
    map<String,DoubleReal>::iterator methyl32 = label_identifiers.find("Methyl32");
    map<String,DoubleReal>::iterator dicpl4 = label_identifiers.find("dICPL4");
    map<String,DoubleReal>::iterator dicpl6 = label_identifiers.find("dICPL6");
    map<String,DoubleReal>::iterator dicpl10 = label_identifiers.find("dICPL10");

    // create string of all labels from advanced section "labels"
    String labels = "Arg6 Arg10 Lys4 Lys6 Lys8 Methyl4 Methyl8 Methyl12 Methyl16 Methyl24 Methyl32 dICPL4 dICPL6 dICPL10";


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


    //--------------------------------------------------
    // section algorithm
    //--------------------------------------------------

    rt_threshold = getParam_().getValue("algorithm:rt_threshold");
    rt_min = getParam_().getValue("algorithm:rt_min");
    intensity_cutoff = getParam_().getValue("algorithm:intensity_cutoff");
    intensity_correlation = getParam_().getValue("algorithm:intensity_correlation");
    model_deviation = getParam_().getValue("algorithm:model_deviation");
    allow_missing_peaks = getFlag_("algorithm:allow_missing_peaks");


    //--------------------------------------------------
    // calculate all possible mass shifts for labelets from section "sample:labels" (concernig missed_cleavage)
    //--------------------------------------------------

    // split string of SILAC labels (selected_labels) and save in a list (SILAClabels)
    vector<String> tempList; // temporary list of strings for SILAC labelets, e.g. "Lys6,Arg8"
    boost::split( tempList, selected_labels, boost::is_any_of("[](){}") ); // any bracket allowed to separate labelets
    for (UInt i = 0; i < tempList.size(); i++)
    {
      if (tempList[i] != "")
      {
        vector<String> tempLabels;
        boost::split( tempLabels, tempList[i], boost::is_any_of(",;: ") ); // various separators allowed to separate labels
        SILAClabels.push_back(tempLabels);
      }
    }

    cout << endl;
    // print SILAC labels
    for (UInt i = 0; i < SILAClabels.size(); i++)
    {
      cout << "SILAC label " << i + 1 << ":   ";
      for (UInt j = 0; j < SILAClabels[i].size(); j++)
      {
        cout << SILAClabels[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;

    // check if all selected labels are included in advanced section "labels"
    for (UInt i = 0; i < SILAClabels.size(); i++)
    {
      for (UInt j = 0; j < SILAClabels[i].size(); ++j)
      {
        Int found = (Int) labels.find(SILAClabels[i][j]);

        if (found < 0)
        {
          throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,SILAClabels[i][j]);
        }
      }
    }

    // generate list of mass shifts
    for (Int ArgPerPeptide = 0; ArgPerPeptide <= missed_cleavages + 1; ArgPerPeptide++)
    {
      for (Int LysPerPeptide = 0; LysPerPeptide <= missed_cleavages + 1; LysPerPeptide++)
      {
        for (Int MethylPerPeptide = 0; MethylPerPeptide <= missed_cleavages + 1; MethylPerPeptide++)
        {
          for (Int dICPLPerPeptide = 0; dICPLPerPeptide <= missed_cleavages + 1; dICPLPerPeptide++)
          {
            if ( ArgPerPeptide + LysPerPeptide + MethylPerPeptide + dICPLPerPeptide > 0 && ArgPerPeptide + LysPerPeptide + MethylPerPeptide + dICPLPerPeptide <= missed_cleavages + 1 )
            {
              vector<DoubleReal> massShiftVector;
              for (UInt i = 0; i < SILAClabels.size(); i++)
              {
                DoubleReal massShift = 0;
                // Considering the case of an amino acid (e.g. LysPerPeptide != 0) for which no label is present (e.g. Lys4There + Lys8There == 0) makes no sense. Therefore each amino acid will have to give its "Go Ahead" before the shift is calculated.
                bool goAhead_Lys = false;
                bool goAhead_Arg = false;
                bool goAhead_Methyl = false;
                bool goAhead_dICPL = false;

                for (UInt j = 0; j < SILAClabels[i].size(); j++)
                {
                  Int Arg6There = 0;	// Is Arg6 in the SILAC label?
                  Int Arg10There = 0;
                  Int Lys4There = 0;
                  Int Lys6There = 0;
                  Int Lys8There = 0;
                  Int Methyl4There = 0;
                  Int Methyl8There = 0;
                  Int Methyl12There = 0;
                  Int Methyl16There = 0;
                  Int Methyl24There = 0;
                  Int Methyl32There = 0;
                  Int dICPL4There = 0;
                  Int dICPL6There = 0;
                  Int dICPL10There = 0;

                  if ( SILAClabels[i][j].find("Arg6") == 0 ) Arg6There = 1;
                  if ( SILAClabels[i][j].find("Arg10") == 0 ) Arg10There = 1;
                  if ( SILAClabels[i][j].find("Lys4") == 0 ) Lys4There = 1;
                  if ( SILAClabels[i][j].find("Lys6") == 0 ) Lys6There = 1;
                  if ( SILAClabels[i][j].find("Lys8") == 0 ) Lys8There = 1;
                  if ( SILAClabels[i][j].find("Methyl4") == 0 ) Methyl4There = 1;
                  if ( SILAClabels[i][j].find("Methyl8") == 0 ) Methyl8There = 1;
                  if ( SILAClabels[i][j].find("Methyl12") == 0 ) Methyl12There = 1;
                  if ( SILAClabels[i][j].find("Methyl16") == 0 ) Methyl16There = 1;
                  if ( SILAClabels[i][j].find("Methyl24") == 0 ) Methyl24There = 1;
                  if ( SILAClabels[i][j].find("Methyl32") == 0 ) Methyl32There = 1;
                  if ( SILAClabels[i][j].find("dICPL4") == 0 ) dICPL4There = 1;
                  if ( SILAClabels[i][j].find("dICPL6") == 0 ) dICPL6There = 1;
                  if ( SILAClabels[i][j].find("dICPL10") == 0 ) dICPL10There = 1;

                  goAhead_Arg = goAhead_Arg || !( (ArgPerPeptide != 0 && Arg6There + Arg10There == 0) );
                  goAhead_Lys = goAhead_Lys || !( (LysPerPeptide != 0 && Lys4There + Lys6There + Lys8There == 0) );
                  goAhead_Methyl = goAhead_Methyl || !( (MethylPerPeptide != 0 && Methyl4There + Methyl8There + Methyl12There + Methyl16There + Methyl24There + Methyl32There == 0) );
                  goAhead_dICPL = goAhead_dICPL || !( (dICPLPerPeptide != 0 && dICPL4There + dICPL6There + dICPL10There == 0) );

                  massShift = massShift + ArgPerPeptide * ( Arg6There * (arg6->second) + Arg10There * (arg10->second) ) + LysPerPeptide * ( Lys4There * (lys4->second) + Lys6There * (lys6->second) + Lys8There * (lys8->second) ) +  MethylPerPeptide * (Methyl4There * (methyl4->second) + Methyl8There * (methyl8->second) + Methyl12There * (methyl12->second) + Methyl16There * (methyl16->second) + Methyl24There * (methyl24->second) + Methyl32There * (methyl32->second) ) + dICPLPerPeptide * ( dICPL4There * (dicpl4->second) + dICPL6There * (dicpl6->second) + dICPL10There * (dicpl10->second) );
                }

                if (goAhead_Arg && goAhead_Lys && goAhead_Methyl && goAhead_dICPL)
                  massShiftVector.push_back(massShift);
              }

              if (!massShiftVector.empty())
                massShifts.push_back(massShiftVector);
            }
          }
        }
      }
    }

    // create zero-mass-shift to search for peptides if no label is specified
    if (massShifts.size() == 0)
    {
      vector<DoubleReal> mass_shift_vector_peptide(1, 0.0);
      massShifts.push_back(mass_shift_vector_peptide);
    }

    // sort the mass shift vector
    sort(massShifts.begin(), massShifts.end());

    // print mass shifts
    for (UInt i = 0; i < massShifts.size(); i++)
    {
      cout << "mass shift " << i + 1 << ":   ";
      for (UInt j = 0; j < massShifts[i].size(); j++)
      {
        cout << massShifts[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }


  //--------------------------------------------------
  // filtering
  //--------------------------------------------------

  void filterData(MSExperiment<Peak1D>& exp, const PeakWidthEstimator::Result &peak_width)
  {
    list<SILACFilter> filters;

    // create filters for all numbers of isotopes per peptide, charge states and mass shifts
    // iterate over all number for peaks per peptide (from max to min)
    for (UInt isotopes_per_peptide = isotopes_per_peptide_max; isotopes_per_peptide >= isotopes_per_peptide_min; isotopes_per_peptide--)
    {
      // iterate over all charge states (from max to min)
      for (UInt charge = charge_max; charge >= charge_min; charge--)
      {
        // iterate over all mass shifts
        for (UInt i = 0; i < massShifts.size(); i++)
        {
          // convert vector<DoubleReal> to set<DoubleReal> for SILACFilter
          vector<DoubleReal> massShifts_set = massShifts[i];

          //copy(massShifts[i].begin(), massShifts[i].end(), inserter(massShifts_set, massShifts_set.end()));
          filters.push_back(SILACFilter(massShifts_set, charge, model_deviation, isotopes_per_peptide, intensity_cutoff, intensity_correlation, allow_missing_peaks));
        }
      }
    }

    // create filtering
    SILACFiltering filtering(exp, peak_width, intensity_cutoff, out_debug);
    filtering.setLogType(log_type_);

    // register filters to the filtering
    for (list<SILACFilter>::iterator filter_it = filters.begin(); filter_it != filters.end(); ++filter_it)
    {
      filtering.addFilter(*filter_it);
    }

    // perform filtering
    filtering.filterDataPoints();

    // retrieve filtered data points
    for (SILACFiltering::Filters::iterator filter_it = filtering.filters_.begin(); filter_it != filtering.filters_.end(); ++filter_it)
    {
      data.push_back(filter_it->getElements());
    }


    //--------------------------------------------------
    // combine DataPoints to improve the clustering
    //--------------------------------------------------

    // DataPoints that originate from filters with same charge state and mass shift(s)
    // and whose filters only differ in number of isotopes per peptide are combined
    // to get one cluster for peptides whose elution profile varies in number of isotopes per peptide

    // perform combination only if the user specified a peaks_per_peptide range > 1
    if (isotopes_per_peptide_min != isotopes_per_peptide_max)
    {
      // erase empty filter results from "data"
      vector<vector<SILACPattern> > data_temp;

      for (vector<vector<SILACPattern> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
      {
        if (data_it->size() != 0)
        {
          data_temp.push_back(*data_it);     // keep DataPoint if it is not empty
        }
      }

      data.swap(data_temp);     // data = data_temp
      data_temp.clear();      // clear "data_temp"

      if (data.size() >= 2)
      {
        Int temp = 0;
        // combine corresponding DataPoints
        vector<vector<SILACPattern> >::iterator data_it_1 = data.begin();      // first iterator over "data" to get first DataPoint for combining
        vector<vector<SILACPattern> >::iterator data_it_2 = data_it_1 + 1;     // second iterator over "data" to get second DataPoint for combining
        vector<vector<SILACPattern> >::iterator data_it_end = data.end() - 1;      // pointer to second last elemnt of "data"
        vector<SILACPattern>::iterator it_1;     // first inner iterator over elements of first DataPoint
        vector<SILACPattern>::iterator it_2;     // second inner iterator over elements of second DataPoint

        while (data_it_1 < data_it_end)      // check for combining as long as first DataPoint is not second last elment of "data"
        {          
          while (data_it_1->size() == 0 && data_it_1 < data_it_end)
          {
            ++data_it_1;      // get next first DataPoint
            data_it_2 = data_it_1 + 1;      // reset second iterator
          }

          if (data_it_1 == data_it_end && data_it_2 == data.end())     // if first iterator points to last element of "data" and second iterator points to end of "data"
          {            
            break;      // stop combining
          }

          while (data_it_2 < data.end() && data_it_2->size() == 0)      // as long as current second DataPoint is empty and second iterator does not point to end of "data"
          {
            ++data_it_2;      // get next second DataPoint
          }

          if (data_it_2 == data.end())      // if second iterator points to end of "data"
          {            
            data_it_2 = data_it_1 + 1;      // reset second iterator
          }

          it_1 = data_it_1->begin();      // set first inner iterator to first element of first DataPoint
          it_2 = data_it_2->begin();      // set second inner iterator to first element of second DataPoint

          // check if DataPoints are not empty
          if (data_it_1->size() != 0 && data_it_2->size() != 0)
          {
            // check if DataPoints have the same charge state and mass shifts
            if (it_1->charge != it_2->charge || it_1->mass_shifts != it_2->mass_shifts)
            {              
              if (data_it_2 < data_it_end)     // if DataPpoints differ and second DataPoint is not second last element of "data"
              {
                temp++;
                ++data_it_2;      // get next second DataPoint
                if (temp > 50000)
                {                  
                  ++data_it_1;
                  temp = 0;
                }
              }

              else if (data_it_2 == data_it_end && data_it_1 < data.end() - 2)     // if DataPpoints differ and second DataPoint is second last element of "data" and first DataPoint is not third last element of "data"
              {
                ++data_it_1;      // get next first DataPoint
                data_it_2 = data_it_1 + 1;      // reset second iterator
              }

              else
              {
                ++data_it_1;      // get next first DataPoint
              }
            }

            else
            {              
              // perform combining
              (*data_it_1).insert(data_it_1->end(), data_it_2->begin(), data_it_2->end());      // append second DataPoint to first DataPoint
              (*data_it_2).clear();     // clear second Datapoint to keep iterators valid and to keep size of "data"

              if (data_it_2 < data_it_end)     // if second DataPoint is not second last element of "data"
              {
                ++data_it_2;      // get next second DataPoint
              }
              else
              {
                data_it_2 = data_it_1 + 1;      // reset second iterator
              }
            }
          }
          else
          {
            ++data_it_1;      // get next first DataPoint
          }
        }

        // erase empty DataPoints from "data"
        vector<vector<SILACPattern> > data_temp;

        for (vector<vector<SILACPattern> >::iterator data_it = data.begin(); data_it != data.end(); ++data_it)
        {
          if (data_it->size() != 0)
          {
            data_temp.push_back(*data_it);     // keep DataPoint if it is not empty
          }
        }

        data.swap(data_temp);     // data = data_temp
        data_temp.clear();      // clear "data_temp"
      }
    }


  }

  ExitCodes main_(int , const char**)
  {
    handleParameters();


    //--------------------------------------------------
    // loading input from .mzML
    //--------------------------------------------------

    MzMLFile file;
    MSExperiment<Peak1D> exp;

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
			std::vector< std::vector< std::pair<String, DoubleReal> > > labels;
			//add none label
			labels.push_back(std::vector< std::pair<String, DoubleReal> > (1,std::make_pair<String,DoubleReal>(String("none"),DoubleReal(0)))); 
			for (Size i = 0; i < SILAClabels.size(); ++i) //SILACLabels MUST be in weight order!!!
			{
				std::vector< std::pair<String, DoubleReal> > one_label;
				for (UInt j = 0; j < SILAClabels[i].size(); ++j)
				{
					one_label.push_back(*(label_identifiers.find( SILAClabels[i][j] ) )); // this dereferencing would break if all SILAClabels would not have been checked before!
				}
				labels.push_back(one_label);
			}
			msq.registerExperiment(exp, labels); //add assays
			msq.assignUIDs();
		}
		MSQuantifications::QUANT_TYPES quant_type = MSQuantifications::MS1LABEL;
		msq.setAnalysisSummaryQuantType(quant_type);//add analysis_summary_

    //--------------------------------------------------
    // estimate peak width
    //--------------------------------------------------

    PeakWidthEstimator::Result peak_width;
    try
    {
      peak_width = estimatePeakWidth(exp);
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

      filterData(exp, peak_width);


      //--------------------------------------------------
      // store filter results
      //--------------------------------------------------

      if (out_filters != "")
      {
        ConsensusMap map;
        for (std::vector<std::vector<SILACPattern> >::const_iterator it = data.begin(); it != data.end(); ++it)
				{        
					generateFilterConsensusByPattern(map, *it);
				}
        writeConsensus(out_filters, map);
      }
    }
    else
    {
      //--------------------------------------------------
      // load filter results
      //--------------------------------------------------

      ConsensusMap map;
      readConsensus(in_filters, map);
      readFilterConsensusByPattern(map);
    }


    //--------------------------------------------------
    // clustering
    //--------------------------------------------------

    clusterData(exp, peak_width);


    //--------------------------------------------------------------
    // write output
    //--------------------------------------------------------------

    if (out_debug != "")
    {
      std::ofstream out((out_debug + ".clusters.csv").c_str());
      
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
        generateClusterDebug(out, **it, cluster_id);
      }
    }
    
    if (out != "")
    {
      ConsensusMap map;
			
      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        generateClusterConsensusByCluster(map, **it);
      }

      // XXX: Need a map per mass shift
      ConsensusMap::FileDescriptions &desc = map.getFileDescriptions();
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


      writeConsensus(out, map);
			if (out_mzq != "")
			{
				ConsensusMap numap(map);
				//calc. ratios
				for (ConsensusMap::iterator cit = numap.begin(); cit != numap.end(); ++cit)
				{
					//~ make ratio templates
					std::vector<ConsensusFeature::Ratio> rts;
					for (std::vector<MSQuantifications::Assay>::const_iterator ait = msq.getAssays().begin()+1; ait != msq.getAssays().end(); ++ait)
					{
						ConsensusFeature::Ratio r;
						r.numerator_ref_ = String(msq.getAssays().begin()->uid_);
						r.denominator_ref_ = String(ait->uid_);				
						r.description_.push_back("Simple ratio calc");
						r.description_.push_back("light to medium/.../heavy");
						//~ "<cvParam cvRef=\"PSI-MS\" accession=\"MS:1001132\" name=\"peptide ratio\"/>"
						rts.push_back(r);
					}
					const std::set< FeatureHandle,FeatureHandle::IndexLess>& feature_handles = cit->getFeatures();
					if (feature_handles.size() > 1)
					{
						std::set< FeatureHandle,FeatureHandle::IndexLess>::const_iterator fit = feature_handles.begin(); // this is unlabeled
						fit++;
						for (fit; fit != feature_handles.end(); ++fit)
						{
							Size ri = std::distance(feature_handles.begin(), fit);
							rts[ri-1].ratio_value_ =  feature_handles.begin()->getIntensity() / fit->getIntensity(); // a proper silacalanyzer algo should never have 0-intensities so no 0devison ...
						}
					}
					
					cit->setRatios(rts);
				}
				msq.addConsensusMap(numap);//add SILACAnalyzer result

				//~ msq.addFeatureMap();//add SILACAnalyzer evidencetrail as soon as clear what is realy contained in the featuremap
				//~ add AuditCollection - no such concept in TOPPTools yet
				writeMzQuantML(out_mzq,msq);
			}
    }

    if (out_clusters != "")
    {
      ConsensusMap map;
      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        UInt cluster_id = 0;
        generateClusterConsensusByPattern(map, **it, cluster_id);
      }

      ConsensusMap::FileDescription &desc = map.getFileDescriptions()[0];
      desc.filename = in;
      desc.label = "Cluster";
      
      writeConsensus(out_clusters, map);
    }

    if (out_features != "")
    {
      FeatureMap<> map;
      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        generateClusterFeatureByCluster(map, **it);
      }

      writeFeatures(out_features, map);
    }

    return EXECUTION_OK;
  }

  void clusterData(const MSExperiment<> &, const PeakWidthEstimator::Result &);

private:
  PeakWidthEstimator::Result estimatePeakWidth(const MSExperiment<Peak1D> &exp);

  /**
   * @brief Generate ConsensusMap from clustering result
   */
  void generateClusterConsensusByCluster(ConsensusMap &, const Clustering &) const;
  
  /**
   * @brief Generate ConsensusMap from clustering result, one consensus per pattern
   */
  void generateClusterConsensusByPattern(ConsensusMap &, const Clustering &, UInt &cluster_id) const;
  
  /**
   * @brief Generate debug output from clustering result
   */
  void generateClusterDebug(std::ostream &out, const Clustering &clustering, UInt &cluster_id) const;
  
  /**
   * @brief Generate ConsensusMap from filter result
   */
  void generateFilterConsensusByPattern(ConsensusMap &, const std::vector<SILACPattern> &) const;
  
  /**
   * @brief Generate a consensus entry from a pattern
   */
  ConsensusFeature generateSingleConsensusByPattern(const SILACPattern &) const;
  
  /**
   * @brief Generate FeatureMap from clustering result
   */
  void generateClusterFeatureByCluster(FeatureMap<> &, const Clustering &) const;
  
  /**
   * @brief Read filter result from ConsensusMap
   */
  void readFilterConsensusByPattern(ConsensusMap &);

  static const String &selectColor(UInt nr);
 
  /**
   * @brief Read consensusXML from file to ConsensusMap
   */
  void readConsensus(const String &filename, ConsensusMap &in) const
  {
    ConsensusXMLFile c_file;
    c_file.load(filename, in);
  }

  /**
   * @brief Write consensusXML from ConsensusMap to file
   */
  void writeConsensus(const String &filename, ConsensusMap &out) const
  {
    out.sortByPosition();
    out.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    out.setExperimentType("silac");

    ConsensusXMLFile c_file;
    c_file.store(filename, out);
  }

	
	  /**
   * @brief Write MzQuantML from ConsensusMap to file
   */
  void writeMzQuantML(const String &filename, MSQuantifications &msq) const
  {
		//~ TODO apply above to ConsensusMap befor putting into Msq
    //~ out.sortByPosition();
    //~ out.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    //~ out.setExperimentType("SILAC");

    MzQuantMLFile file;
    file.store(filename, msq);
  }
	
  /**
   * @brief Write featureXML from FeatureMap to file
   */
  void writeFeatures(const String &filename, FeatureMap<> &out) const
  {
    out.sortByPosition();
    out.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    FeatureXMLFile f_file;
    f_file.store(filename, out);
  }
};

void TOPPSILACAnalyzer::clusterData(const MSExperiment<> &exp, const PeakWidthEstimator::Result &peak_width)
{
  typedef Clustering::PointCoordinate PointCoordinate;

  ProgressLogger progresslogger;
  progresslogger.setLogType(log_type_);
  progresslogger.startProgress(0, data.size(), "clustering data");

  // Use peak half width @1000 Th for mz threshold
  DoubleReal mz_threshold = peak_width(1000);

  // Use double median of spectrum spacing for max rt spacing
  DoubleReal rt_max_spacing = 0;
  {
    // Calculate distance between each spectrum; this needs sorted spectra
    std::vector<Real> space;
    MSExperiment<>::const_iterator it1 = exp.begin();
    MSExperiment<>::const_iterator it2 = exp.begin(); ++it2;
    for (; it2 != exp.end(); ++it1, ++it2)
    {
      DoubleReal s = it2->getRT() - it1->getRT();
      space.push_back(s);
    }

    sort(space.begin(), space.end());

    // Calculate median by extracting the middle element (okay, the upper median)
    // Set max spacing to five times the median spectrum spacing
    // The five is an empirical value
    if (space.size()) rt_max_spacing = space[space.size() / 2 + 1] * 5;
  }

  UInt data_id = 0;

  for (vector<vector<SILACPattern> >::iterator data_it = data.begin();
       data_it != data.end();
       ++data_it, ++data_id)
  {
    const PointCoordinate max_delta(rt_threshold, mz_threshold);
    Clustering *clustering = new Clustering(max_delta, rt_min, rt_max_spacing);

    for (vector<SILACPattern>::iterator it = data_it->begin(); it != data_it->end(); ++it)
    {
      const PointCoordinate key(it->rt, it->mz);
      SILACPattern &p = *it;
      clustering->insertPoint(key, &p);
    }

    clustering->cluster();

    cluster_data.push_back(clustering);

    progresslogger.setProgress(data_id);
  }

  progresslogger.endProgress();
}

PeakWidthEstimator::Result TOPPSILACAnalyzer::estimatePeakWidth(const MSExperiment<Peak1D> &exp)
{
  ProgressLogger progresslogger;
  progresslogger.setLogType(log_type_);
  progresslogger.startProgress(0, 1, "estimate peak width");

  PeakWidthEstimator::Result ret = PeakWidthEstimator::estimateFWHM(exp);

  progresslogger.endProgress();
  std::cout << "Estimated peak width: e ^ (" << ret.c0 << " + " << ret.c1 << " * log mz)" << std::endl;
  return ret;
}

void TOPPSILACAnalyzer::generateClusterConsensusByCluster(ConsensusMap &out, const Clustering &clustering) const
{
  // iterate over clusters
  for (Clustering::Grid::const_iterator cluster_it = clustering.grid.begin(); cluster_it != clustering.grid.end(); ++cluster_it)
  {
    ConsensusFeature consensus;

    // determine the number of peptides
    // for that we look at the first point in the first pattern
    const SILACPattern &firstPattern = *(cluster_it->second.begin())->second;
    const SILACPoint &firstPoint = *(firstPattern.points.begin());
    UInt numberPeptides = firstPoint.intensities.size();
    UInt charge = firstPoint.charge;

    // sums for each peptide of the pair (triplet, singlet, ...)
    std::vector<DoubleReal> sumMzIntensities (numberPeptides,0);    // sum m/z * intensity (for intensity-weighted m/z average)
    std::vector<DoubleReal> sumRtIntensities (numberPeptides,0);    // sum rt * intensity (for intensity-weighted rt average)
    std::vector<DoubleReal> sumIntensities (numberPeptides,0);    // sum intensity (for 'feature volume' = peptide intensity)
    std::vector<DoubleReal> sumIntensitiesMonoisotopic (numberPeptides,0);    // sum intensity of monoisotopic mass trace (for normalisation of intensity-weighted m/z average)
    std::vector<DoubleReal> maxIntensityXIC (numberPeptides,0);    // tracks maximum of sumIntensitiesXIC
    std::vector<DoubleReal> RtAtMaxIntensityXIC (numberPeptides,0);    // tracks rt at maximum of sumIntensitiesXIC
    
    // iterate over SILAC patterns in each cluster
    for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin();
         pattern_it != cluster_it->second.end();
         ++pattern_it)
    {
      const SILACPattern &pattern = *pattern_it->second;
      std::vector<DoubleReal> sumIntensitiesXIC (numberPeptides,0);    // sums intensities at fixed rt (for XIC)

      // iterate over SILAC points in each SILAC pattern
      for (std::vector<SILACPoint>::const_iterator point_it = pattern.points.begin();
           point_it != pattern.points.end();
           ++point_it)
      {
        const SILACPoint &point = *point_it;

        // iterate over peptides in doublet (or triplet, ...)
        UInt peptide = 0;    // peptide for which intensity, retention time and m/z are to be calculated
        for (std::vector<std::vector<DoubleReal> >::const_iterator peptide_it = point.intensities.begin();
             peptide_it != point.intensities.end();
             ++peptide_it)
      {
          // iterate over isotopes in peptide
          UInt isotope = 0;
          for (std::vector<DoubleReal>::const_iterator isotope_it = peptide_it->begin();
               isotope_it != peptide_it->end();
               ++isotope_it)
        {
            sumIntensities[peptide] += point.intensities[peptide][isotope];
            sumIntensitiesXIC[peptide] += point.intensities[peptide][isotope];
            sumRtIntensities[peptide] += point.rt * point.intensities[peptide][isotope];
            if (isotope == 0)
            {
              sumMzIntensities[peptide] += pattern.mz_positions[peptide][isotope] * point.intensities[peptide][isotope];
              sumIntensitiesMonoisotopic[peptide] += point.intensities[peptide][isotope];
            }
            //sumMzIntensities[peptide] += (pattern.mz_positions[peptide][isotope] - (isotope * 1.003355 / point.charge)) * point.intensities[peptide][isotope];
            ++isotope;
          }
          ++peptide;
        }

        }
      
      // check for each peptide if its XIC intensity has been raised
      for (UInt peptide = 0; peptide < numberPeptides; ++peptide)
      {
        if (sumIntensitiesXIC[peptide] > maxIntensityXIC[peptide])
        {
          maxIntensityXIC[peptide] = sumIntensitiesXIC[peptide];
          RtAtMaxIntensityXIC[peptide] = pattern.rt;
      }
    }
    }
    /*cout << "light m/z: " << sumMzIntensities[0]/sumIntensitiesMonoisotopic[0] << '\n';
     cout << "light rt (intensity averaged): " << sumRtIntensities[0]/sumIntensities[0] << '\n';
     cout << "light rt (at max XIC): " << RtAtMaxIntensityXIC[0] << '\n';
     cout << "heavy m/z: " << sumMzIntensities[1]/sumIntensitiesMonoisotopic[1] << '\n';
     cout << "heavy rt (intensity averaged): " << sumRtIntensities[1]/sumIntensities[1] << '\n';
     cout << "heavy rt (at max XIC): " << RtAtMaxIntensityXIC[1] << '\n' << '\n';*/

    // consensus feature has coordinates of the light peptide
    consensus.setMZ(sumMzIntensities[0]/sumIntensitiesMonoisotopic[0]);    // intensity-average only over the mono-isotopic peak
    consensus.setRT(sumRtIntensities[0]/sumIntensities[0]);    // intensity-average over the entire peptide, i.e. all peptides
    //consensus.setRT(RtAtMaxIntensityXIC[0]);
    consensus.setIntensity(sumIntensities[0]);
    consensus.setCharge(charge);
    consensus.setQuality(std::floor(firstPattern.mass_shifts[1] * charge));    // set Quality to the first mass shift (allows later to filter in consensXML)

    // attach features to consensus
    for (UInt peptide = 0; peptide < numberPeptides; ++peptide)
    {
      FeatureHandle feature;

      feature.setMZ(sumMzIntensities[peptide]/sumIntensitiesMonoisotopic[peptide]);
      feature.setRT(sumRtIntensities[peptide]/sumIntensities[peptide]);
      //feature.setRT(RtAtMaxIntensityXIC[peptide]);
      feature.setIntensity(sumIntensities[peptide]);
      feature.setCharge(charge);
      feature.setMapIndex(peptide);
      out.getFileDescriptions()[peptide].size++;

      consensus.insert(feature);
    }

    // add consensus to consensus map
    out.push_back(consensus);
  }  
}

void TOPPSILACAnalyzer::generateClusterConsensusByPattern(ConsensusMap &out, const Clustering &clustering, UInt &cluster_id) const
    {
  for (Clustering::Grid::const_iterator cluster_it = clustering.grid.begin(); cluster_it != clustering.grid.end(); ++cluster_it, ++cluster_id)
  {
    for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin(); pattern_it != cluster_it->second.end(); ++pattern_it)
    {
      ConsensusFeature consensus = generateSingleConsensusByPattern(*pattern_it->second);

      consensus.setMetaValue("color", selectColor(cluster_id));
      consensus.setMetaValue("Cluster ID", cluster_id);

      out.getFileDescriptions()[0].size++;

      out.push_back(consensus);
    }
  }
}

void TOPPSILACAnalyzer::generateClusterDebug(std::ostream &out, const Clustering &clustering, UInt &cluster_id) const
{
  for (Clustering::Grid::const_iterator cluster_it = clustering.grid.begin();
       cluster_it != clustering.grid.end();
       ++cluster_it, ++cluster_id)
  {
      for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin();
           pattern_it != cluster_it->second.end();
           ++pattern_it)
      {
      const SILACPattern &pattern = *pattern_it->second;

      std::ostringstream preamble;

      preamble
      << std::fixed << std::setprecision(8)
      << cluster_id << ','
      << pattern.rt << ','
      << pattern.mz << ','
      << pattern.charge << ',';

      for (std::vector<DoubleReal>::const_iterator shift_it = pattern.mass_shifts.begin();
           shift_it != pattern.mass_shifts.end();
           ++shift_it)
        {
        preamble
        << *++shift_it * pattern.charge << ',';
        }
      
      for (std::vector<std::vector<DoubleReal> >::const_iterator shift_inten_it = pattern.intensities.begin();
           shift_inten_it != pattern.intensities.end();
           ++shift_inten_it)
      {
        UInt peak_inten_id = 0;
        for (std::vector<DoubleReal>::const_iterator peak_inten_it = shift_inten_it->begin();
             peak_inten_it != shift_inten_it->end();
             ++peak_inten_it, ++peak_inten_id)
        {
          preamble
          << *peak_inten_it << ',';
      }
        for (; peak_inten_id < isotopes_per_peptide_max; ++peak_inten_id)
        {
          preamble
          << "NA,";
        }
      }

      for (std::vector<SILACPoint>::const_iterator point_it = pattern.points.begin();
           point_it != pattern.points.end();
           ++point_it)
      {
        const SILACPoint &point = *point_it;

        out
        << preamble.str()
        << point.mz;

        // write INT_RAW_...
        for (std::vector<std::vector<DoubleReal> >::const_iterator shift_inten_it = point.intensities.begin();
             shift_inten_it != point.intensities.end();
             ++shift_inten_it)
        {
          UInt peak_inten_id = 0;
          for (std::vector<DoubleReal>::const_iterator peak_inten_it = shift_inten_it->begin();
               peak_inten_it != shift_inten_it->end();
               ++peak_inten_it, ++peak_inten_id)
          {
            out << ','
            << *peak_inten_it;
    }
          for (; peak_inten_id < isotopes_per_peptide_max; ++peak_inten_id)
          {
            out
            << ",NA";
  }
}

        // write MZ_RAW_...
        for (std::vector<std::vector<DoubleReal> >::const_iterator shift_mz_it = point.mz_positions.begin();
             shift_mz_it != point.mz_positions.end();
             ++shift_mz_it)
{
          UInt peak_mz_id = 0;
          for (std::vector<DoubleReal>::const_iterator peak_mz_it = shift_mz_it->begin();
               peak_mz_it != shift_mz_it->end();
               ++peak_mz_it, ++peak_mz_id)
  {
            out << ','
            << *peak_mz_it;
          }
          for (; peak_mz_id < isotopes_per_peptide_max; ++peak_mz_id)
    {
            out
            << ",NA";
          }
        }

        out << '\n';
    }
  }
}
}

void TOPPSILACAnalyzer::generateFilterConsensusByPattern(ConsensusMap &out, const std::vector<SILACPattern> &pattern) const
{
  for (std::vector<SILACPattern>::const_iterator pattern_it = pattern.begin(); pattern_it != pattern.end(); ++pattern_it)
  {
    out.push_back(generateSingleConsensusByPattern(*pattern_it));
  }
}

ConsensusFeature TOPPSILACAnalyzer::generateSingleConsensusByPattern(const SILACPattern &pattern) const
{
  // XXX: get from experiment
  Int charge = pattern.charge;

  ConsensusFeature consensus;
  consensus.setRT(pattern.rt);
  consensus.setMZ(pattern.mz);
  consensus.setIntensity(pattern.intensities[0][0]);
  consensus.setCharge(charge);

  consensus.setMetaValue("Peaks per peptide", pattern.isotopes_per_peptide);

  // Output mass shifts
  {
    std::ostringstream out;
    out << std::fixed << std::setprecision(4);
    for (vector<DoubleReal>::const_iterator shift_it = pattern.mass_shifts.begin() + 1; shift_it != pattern.mass_shifts.end(); ++shift_it)
    {
      out << *shift_it * charge << ';';
    }
    // Remove the last delimiter
    std::string outs = out.str(); outs.erase(outs.end() - 1);
    consensus.setQuality(std::floor(pattern.mass_shifts.at(1) * charge));
    consensus.setMetaValue("SILAC", outs);
  }

  // Output all intensities per peptide as list
  {
    std::ostringstream out;
    for (vector<vector<DoubleReal> >::const_iterator inten_it = pattern.intensities.begin(); inten_it != pattern.intensities.end(); ++inten_it)
    {
      std::ostringstream out2;
      out2 << std::fixed << std::setprecision(4);
      for (vector<DoubleReal>::const_iterator inten2_it = inten_it->begin(); inten2_it != inten_it->end(); ++inten2_it)
      {
        out2 << *inten2_it << ',';
      }
      // Remove the last delimiter
      std::string out2s = out2.str(); out2s.erase(out2s.end() - 1);
      out << out2s << ';';
    }
    // Remove the last delimiter
    std::string outs = out.str(); outs.erase(outs.end() - 1);
    consensus.setMetaValue("Intensities", outs);
  }

  UInt point_id = 0;
  for (std::vector<SILACPoint>::const_iterator point_it = pattern.points.begin();
       point_it != pattern.points.end();
       ++point_it, ++point_id)
  {
    FeatureHandle point;
    point.setRT(point_it->rt);
    point.setMZ(point_it->mz);
    point.setUniqueId(point_id);

    consensus.insert(point);
  }

  return consensus;
}

void TOPPSILACAnalyzer::generateClusterFeatureByCluster(FeatureMap<> &out, const Clustering &clustering) const
{
  for (Clustering::Grid::const_iterator cluster_it = clustering.grid.begin(); cluster_it != clustering.grid.end(); ++cluster_it)
  {
    // RT value as weighted RT position of all peaks
    DoubleReal global_rt = 0;
    // Total intensity
    DoubleReal global_intensity = 0;

    for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin();
         pattern_it != cluster_it->second.end();
         ++pattern_it)
    {
      SILACPattern &pattern = *pattern_it->second;

      for (std::vector<std::vector<DoubleReal> >::const_iterator shift_inten_it = pattern.intensities.begin();
           shift_inten_it != pattern.intensities.end();
           ++shift_inten_it)
      {
        for (std::vector<DoubleReal>::const_iterator peak_inten_it = shift_inten_it->begin();
             peak_inten_it != shift_inten_it->end();
             ++peak_inten_it)
        {
          DoubleReal intensity = *peak_inten_it;

          // Add to RT value and global intensity
          global_rt += intensity * pattern.rt;
          global_intensity += intensity;
        }
      }
    }

    // Calculate global RT value
    global_rt /= global_intensity;

    SILACPattern &pattern_first = *cluster_it->second.begin()->second;

    for (UInt shift_id = 0; shift_id < pattern_first.mass_shifts.size(); ++shift_id)
    {
      // XXX: Feature detection produces a stray 0 mass shift
      if (shift_id > 0 && pattern_first.mass_shifts[shift_id] == 0)
        continue;

      Feature feature;

      // MZ value as weighted MZ position of monoisotopic peaks of given mass shift
      DoubleReal shift_mz = 0;
      // Total intensity
      DoubleReal shift_intensity = 0;
      // Total intensity of monoisotopic peak
      DoubleReal shift_intensity0 = 0;

      // Bounding box per peak
      std::map<UInt, DBoundingBox<2> > bboxs;

      for (Clustering::Cluster::const_iterator pattern_it = cluster_it->second.begin();
           pattern_it != cluster_it->second.end();
           ++pattern_it)
      {
        SILACPattern &pattern = *pattern_it->second;

        const std::vector<DoubleReal> &intensities = pattern.intensities[shift_id];
        DoubleReal mz = pattern.mz + pattern.mass_shifts[shift_id];
        DoubleReal intensity0 = intensities[0];

        // Add to MZ value and shift intensity of monoisotopic peak
        shift_mz += intensity0 * mz;
        shift_intensity0 += intensity0;

        // Iterator over every peak
        UInt peak_id = 0;
        std::vector<DoubleReal>::const_iterator peak_inten_it = intensities.begin();
        DoubleReal peak_mz = mz;
        for (;
             peak_inten_it != intensities.end();
             ++peak_id, ++peak_inten_it, peak_mz += 1. / pattern.charge)
        {
          shift_intensity += *peak_inten_it;
          bboxs[peak_id].enlarge(pattern.rt, peak_mz);
        }
      }

      // Add each bbox as convex hulls to the cluster
      for (std::map<UInt, DBoundingBox<2> >::const_iterator bboxs_it = bboxs.begin();
           bboxs_it != bboxs.end();
           ++bboxs_it)
      {
        ConvexHull2D hull;
        hull.addPoint(bboxs_it->second.min_);
        hull.addPoint(bboxs_it->second.max_);
        feature.getConvexHulls().push_back(hull);
      }

      // XXX: Real quality?
      feature.setOverallQuality(1);
      feature.setCharge(pattern_first.charge);

      // Calculate MZ value
      shift_mz /= shift_intensity0;

      feature.setRT(global_rt);
      feature.setMZ(shift_mz);
      feature.setIntensity(shift_intensity);

      out.push_back(feature);
    }
  }
}

void TOPPSILACAnalyzer::readFilterConsensusByPattern(ConsensusMap &in)
{
  std::map<std::pair<Int, Int>, std::vector<SILACPattern> > layers;

  for (ConsensusMap::const_iterator pattern_it = in.begin(); pattern_it != in.end(); ++pattern_it)
  {
    SILACPattern pattern;
    pattern.rt = pattern_it->getRT();
    pattern.mz = pattern_it->getMZ();
    pattern.charge = pattern_it->getCharge();
    pattern.quality = pattern_it->getQuality();

    pattern.isotopes_per_peptide = pattern_it->getMetaValue("Peaks per peptide");

    StringList text = StringList::create(pattern_it->getMetaValue("Mass shifts [Da]"), ';');
    pattern.mass_shifts.push_back(0);
    for (StringList::const_iterator text_it = text.begin(); text_it != text.end(); ++text_it)
    {
      pattern.mass_shifts.push_back(text_it->toDouble() / pattern.charge);
    }

    text = StringList::create(pattern_it->getMetaValue("Intensities"), ';');
    for (StringList::const_iterator text_it = text.begin(); text_it != text.end(); ++text_it)
    {
      StringList text2 = StringList::create(*text_it, ',');
      vector<DoubleReal> inten;
      for (StringList::const_iterator text2_it = text2.begin(); text2_it != text2.end(); ++text2_it)
      {
        inten.push_back(text2_it->toDouble());
      }
      pattern.intensities.push_back(inten);
    }

    for (ConsensusFeature::const_iterator point_it = pattern_it->begin(); point_it != pattern_it->end(); ++point_it)
    {
      SILACPoint point;
      point.rt = point_it->getRT();
      point.mz = point_it->getMZ();

      pattern.points.push_back(point);
    }

    layers[std::make_pair(Int(pattern.mass_shifts.at(1)), pattern.charge)].push_back(pattern);
  }

  for (std::map<std::pair<Int, Int>, std::vector<SILACPattern> >::iterator it = layers.begin(); it != layers.end(); ++it)
  {
    data.push_back(it->second);
  }
}

const String &TOPPSILACAnalyzer::selectColor(UInt nr)
{
  // 15 HTML colors
  const static String colors[] = {
    "#00FFFF", "#000000", "#0000FF", "#FF00FF", "#008000",
    "#808080", "#00FF00", "#800000", "#000080", "#808000",
    "#800080", "#FF0000", "#C0C0C0", "#008080", "#FFFF00",
  };
  const Int colors_len = 15;

  return colors[nr % colors_len];
}

//@endcond

int main(int argc, const char** argv )
{
  TOPPSILACAnalyzer tool;
  return tool.main(argc, argv);
}
