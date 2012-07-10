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
  @page TOPP_FeatureFinderRaw FeatureFinderRaw

  @brief Identifies peptide features in raw (i.e. profile) LC-MS data.

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderRaw \f$ \longrightarrow \f$</td>
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

  FeatureFinderRaw is a tool for the identification of peptide features in profile LC-MS data.

  <b>Algorithm</b>

  The underlying algorithm of the tool is equivalent to that of SILACAnalyzer.

  <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_FeatureFinderRaw.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_FeatureFinderRaw.html

  <b>Parameter Tuning</b>

  <i>input:</i>
  - in [*.mzML] - LC-MS dataset to be analyzed
  - ini [*.ini] - file containing all parameters (see discussion below)

  <i>standard output:</i>
 - out [*.consensusXML] - contains the list of identified peptides

 <i>optional output:</i>
 - out_clusters [*.consensusXML] - contains the complete set of data points passing the filters
 
  The results of an analysis can easily visualized within TOPPView. Simply load *.consensusXML and *.featureXML as layers over the original *.mzML.

  Parameters in section <i>algorithm:</i>
 - <i>allow_missing_peaks</i> - Low intensity peaks might be missing from the isotopic pattern of some of the peptides. Specify if such peptides should be included in the analysis.
  - <i>rt_threshold</i> - Upper bound for the retention time [s] over which a characteristic peptide elutes.
  - <i>rt_min</i> - Lower bound for the retentions time [s].
  - <i>intensity_cutoff</i> - Lower bound for the intensity of isotopic peaks in a SILAC pattern.
  - <i>intensity_correlation</i> - Lower bound for the Pearson correlation coefficient, which measures how well intensity profiles of different isotopic peaks correlate.
  - <i>model_deviation</i> - Upper bound on the factor by which the ratios of observed isotopic peaks are allowed to differ from the ratios of the theoretic averagine model, i.e. ( theoretic_ratio / model_deviation ) < observed_ratio < ( theoretic_ratio * model_deviation ).

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderRaw
: public TOPPBase
{
  private:

    // input and output files
    String in;
    String out;

    // section "sample"
    Int charge_min;
    Int charge_max;
    Int missed_cleavages;
    Int isotopes_per_peptide_min;
    Int isotopes_per_peptide_max;

    // section "algorithm"
    DoubleReal rt_threshold;
    DoubleReal rt_min;
    DoubleReal intensity_cutoff;
    DoubleReal intensity_correlation;
    DoubleReal model_deviation;

    vector<vector <DoubleReal> > massShifts;      // list of mass shifts

    typedef SILACClustering Clustering;

    vector<vector<SILACPattern> > data;
    vector<Clustering *> cluster_data;


  public:
    TOPPFeatureFinderRaw()
      : TOPPBase("FeatureFinderRaw","Determination of peak ratios in LC-MS data", true)
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
    // create flag for output file (.featureXML)
    registerOutputFile_("out", "<file>", "", "Set of all identified peptides. The m/z-RT positions correspond to the lightest peptide in each group.", false);
    setValidFormats_("out", StringList::create("featureXML"));

    // create section "sample" for adjusting sample parameters
    registerSubsection_("sample", "Parameters describing the sample and its labels.");
    // create section "algorithm" for adjusting algorithm parameters
    registerSubsection_("algorithm", "Parameters for the algorithm.");
  }


  // create prameters for sections (set default values and restrictions)
  Param getSubsectionDefaults_(const String& section) const
  {
    Param defaults;


    //--------------------------------------------------
    // section sample
    //--------------------------------------------------

    if (section == "sample")
    {
      defaults.setValue("charge", "2:4", "Range of charge states in the sample, i.e. min charge : max charge.");
      defaults.setValue("peaks_per_peptide", "3:5", "Range of peaks per peptide in the sample, i.e. min peaks per peptide : max peaks per peptide. For example 3:6, if isotopic peptide patterns in the sample consist of either three, four, five or six isotopic peaks. ", StringList::create("advanced"));
    }


    //--------------------------------------------------
    // section algorithm
    //--------------------------------------------------

    if (section == "algorithm")
    {
      defaults.setValue("rt_threshold", 30.0, "Typical retention time [s] over which a characteristic peptide elutes. (This is not an upper bound. Peptides that elute for longer will be reported.)");
      defaults.setMinFloat("rt_threshold", 0.0);
      defaults.setValue("rt_min", 0.0, "Lower bound for the retention time [s].", StringList::create("advanced"));
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
    // get name of output file (.featureXML)
    out = getStringOption_("out");


    //--------------------------------------------------
    // section sample
    //--------------------------------------------------

    // get selected missed_cleavages
    missed_cleavages = 0;

    // get selected charge range
    String charge_string = getParam_().getValue("sample:charge");
    DoubleReal charge_min_temp, charge_max_temp;
    parseRange_(charge_string, charge_min_temp, charge_max_temp);
    charge_min = (Int)charge_min_temp;
    charge_max = (Int)charge_max_temp;

    // check if charge_min is smaller than charge max, if not swap
    if (charge_min > charge_max)
      swap(charge_min, charge_max);

    // get selected peaks range
    String isotopes_per_peptide_string = getParam_().getValue("sample:peaks_per_peptide");
    DoubleReal isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp;
    parseRange_(isotopes_per_peptide_string, isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp);
    isotopes_per_peptide_min = (Int)isotopes_per_peptide_min_temp;
    isotopes_per_peptide_max = (Int)isotopes_per_peptide_max_temp;

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


    {
      vector<DoubleReal> mass_shift_vector_peptide(1, 0.0);
      massShifts.push_back(mass_shift_vector_peptide);
    }
  }


  //--------------------------------------------------
  // filtering
  //--------------------------------------------------

  void filterData(MSExperiment<Peak1D>& exp, const PeakWidthEstimator::Result &peak_width)
  {
    list<SILACFilter> filters;

    // create filters for all numbers of isotopes per peptide, charge states and mass shifts
    // iterate over all number for peaks per peptide (from max to min)
    for (Int isotopes_per_peptide = isotopes_per_peptide_max; isotopes_per_peptide >= isotopes_per_peptide_min; isotopes_per_peptide--)
    {
      // iterate over all charge states (from max to min)
      for (Int charge = charge_max; charge >= charge_min; charge--)
      {
        // iterate over all mass shifts
        for (UInt i = 0; i < massShifts.size(); i++)
        {
          // convert vector<DoubleReal> to set<DoubleReal> for SILACFilter
          vector<DoubleReal> massShifts_set = massShifts[i];

          //copy(massShifts[i].begin(), massShifts[i].end(), inserter(massShifts_set, massShifts_set.end()));
          filters.push_back(SILACFilter(massShifts_set, charge, model_deviation, isotopes_per_peptide, intensity_cutoff, intensity_correlation, 0));
        }
      }
    }

    // create filtering
    SILACFiltering filtering(exp, peak_width, intensity_cutoff, "");
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
        if ( !data_it->empty() )
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
          while ( data_it_1->empty() && data_it_1 < data_it_end )
          {
            ++data_it_1;      // get next first DataPoint
            data_it_2 = data_it_1 + 1;      // reset second iterator
          }

          if (data_it_1 == data_it_end && data_it_2 == data.end())     // if first iterator points to last element of "data" and second iterator points to end of "data"
          {            
            break;      // stop combining
          }

          while ( data_it_2 < data.end() && data_it_2->empty() )      // as long as current second DataPoint is empty and second iterator does not point to end of "data"
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
          if ( !data_it_1->empty() && !data_it_2->empty() )
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
          if ( !data_it->empty() )
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


    //--------------------------------------------------
    // filter input data
    //--------------------------------------------------

    filterData(exp, peak_width);


    //--------------------------------------------------
    // clustering
    //--------------------------------------------------

    clusterData(peak_width);


    //--------------------------------------------------------------
    // write output
    //--------------------------------------------------------------

    if (out != "")
    {
      FeatureMap<> map;
      for (vector<Clustering *>::const_iterator it = cluster_data.begin(); it != cluster_data.end(); ++it)
      {
        generateClusterFeatureByCluster(map, **it);
      }

      writeFeatures(out, map);
    }

    return EXECUTION_OK;
  }

  void clusterData(const PeakWidthEstimator::Result &);

private:
  PeakWidthEstimator::Result estimatePeakWidth(const MSExperiment<Peak1D> &);

  void generateClusterFeatureByCluster(FeatureMap<> &, const Clustering &) const;

  void writeFeatures(const String &filename, FeatureMap<> &out) const
  {
    out.sortByPosition();
    out.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    FeatureXMLFile f_file;
    f_file.store(filename, out);
  }
};

void TOPPFeatureFinderRaw::clusterData(const PeakWidthEstimator::Result &peak_width)
{
  typedef Clustering::PointCoordinate PointCoordinate;

  ProgressLogger progresslogger;
  progresslogger.setLogType(log_type_);
  progresslogger.startProgress(0, data.size(), "clustering data");

  // Use peak half width @1000 Th for mz threshold
  DoubleReal mz_threshold = peak_width(1000);

  UInt data_id = 0;

  for (vector<vector<SILACPattern> >::iterator data_it = data.begin();
       data_it != data.end();
       ++data_it, ++data_id)
  {
    const PointCoordinate max_delta(rt_threshold, mz_threshold);
    Clustering *clustering = new Clustering(max_delta, rt_min, 0);

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

PeakWidthEstimator::Result TOPPFeatureFinderRaw::estimatePeakWidth(const MSExperiment<Peak1D> &exp)
{
  ProgressLogger progresslogger;
  progresslogger.setLogType(log_type_);
  progresslogger.startProgress(0, 1, "estimate peak width");

  PeakWidthEstimator::Result ret = PeakWidthEstimator::estimateFWHM(exp);

  progresslogger.endProgress();
  std::cout << "Estimated peak width: e ^ (" << ret.c0 << " + " << ret.c1 << " * log mz)" << std::endl;
  return ret;
}

void TOPPFeatureFinderRaw::generateClusterFeatureByCluster(FeatureMap<> &out, const Clustering &clustering) const
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

//@endcond

int main(int argc, const char** argv )
{
  TOPPFeatureFinderRaw tool;
  return tool.main(argc, argv);
}
