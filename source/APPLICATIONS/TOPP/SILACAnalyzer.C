// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

//OpenMS includes
#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

//clustering
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/COMPARISON/CLUSTERING/AverageLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/CompleteLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
//~ #include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>

//Contrib includes
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_SILACAnalyzer SILACAnalyzer

  @brief Identifies peptide pairs in LC-MS data and determines their relative abundance.

  SILACAnalyzer is a tool for the fully automated analysis of quantitative proteomics data. It identifies pairs of isotopic envelopes with fixed m/z separation. It requires no prior sequence identification of the peptides. In what follows we first explain the algorithm and then discuss the tuning of its parameters.

  <b>Algorithm</b>

  The algorithm is divided into three parts: filtering, clustering and linear fitting, see Fig. (d), (e) and (f). In the following discussion let us consider a particular mass spectrum at retention time 1350 s, see Fig. (a). It contains a peptide of mass 1492 Da and its 6 Da heavier labelled counterpart. Both are doubly charged in this instance. Their isotopic envelopes therefore appear at 746 and 749 in the spectrum. The isotopic peaks within each envelope are separated by 0.5. The spectrum was recorded at finite intervals. In order to read accurate intensities at arbitrary m/z we spline-fit over the data, see Fig. (b).

  We would like to search for such peptide pairs in our LC-MS data set. As a warm-up let us consider a standard intensity cut-off filter, see Fig. (c). Scanning through the entire m/z range (red dot) only data points with intensities above a certain threshold pass the filter. Unlike such a local filter, the filter used in our algorithm takes intensities at a range of m/z positions into account, see Fig. (d). A data point (red dot) passes if
  - all six intensities at m/z, m/z+0.5, m/z+1, m/z+3, m/z+3.5 and m/z+4 lie above a certain threshold and
  - the intensities within the first envelope (at m/z, m/z+0.5 and m/z+1) and second envelope (at m/z+3, m/z+3.5 and m/z+4) decrease successively.

  Let us now filter not only a single spectrum but all spectra in our data set. Data points that pass the filter form clusters in the t-m/z plane, see Fig. (e). Each cluster centers around the unlabelled peptide of a pair. We now use hierarchical clustering methods to assign each data point to a specific cluster. The optimum number of clusters is determined by maximizing the silhouette width of the partitioning. Each data point in a cluster corresponds to three pairs of intensities (at [m/z, m/z+3], [m/z+0.5, m/z+3.5] and [m/z+1, m/z+4]). A plot of all intensity pairs in a cluster shows a clear linear correlation, see Fig. (f). Using linear regression we can determine the relative amounts of labelled and unlabelled peptides in the sample.

  @image html SILACAnalyzer_algorithm.png

  <b>Parameter Tuning</b>
 
  SILACAnalyzer can either search for SILAC pairs (as described in the above paragraph) or triplets:
  - type - double (for SILAC pairs), triple (for SILAC triplets)

  <i>input:</i>
  - in [*.mzML] - LC-MS dataset to be analyzed
  - ini [*.ini] - file containing all parameters (see discussion below)

  <i>standard output:</i>
  - out [*.consensusXML] - contains the list of identified peptide pairs (retention time and m/z of the lighter peptide, heavy-to-light ratio)
  - out_visual [*.featureXML] - contains the complete set of data points (retention time, m/z, intensity) of all peptide pairs

  The results of an analysis can easily visualized within TOPPView. Simply load *.consensusXML and *.featureXML as layers over the original *.mzML.

  <i>optional output:</i>
  @n If -silac_debug is enabled, SILACAnalyzer generates a number of further files:
  - [*.dat] - contains the list of identified peptide pairs in a simple text file, c.f. *.consensusXML
  - [*_clusters.dat] -  contains the complete set of data points of all peptide pairs in a simple text file, c.f. *.featureXML
  - [*.input] - gnuplot script for the visualization of the results. Running (gnuplot *.input) generates a number of *.eps plots. The range of clusters to be plotted can be specified by the parameters -cluster_min/max.

  The following parameters are straightforward:
  - mass_separation_light_medium - mass gap between light and medium isotopic envelopes [Da] (only relevant for the search for SILAC triplets, i.e. type triple)
  - mass_separation_light_heavy - mass gap between light and heavy isotopic envelopes [Da]
  - charge_min/max - range of charge states
  - mz_step_width - step width with which the interpolated spectrum, Fig. (b), is scanned. The step width should be of about the same order with which the raw data were recorded, see Fig. (a).

  The remaining parameters should be tuned in the following order:
  - intensity_cutoff - adjust the intensity cutoff such that the data points that pass the non-local filter (*.featureXML layer) form clear distinct clusters,  see Fig. (e). Ignore the coloring of the clusters at that stage.
  - rt_scaling - pick a representative cluster. rt_scaling = (width of the cluster in Da)/(height of the cluster in sec)
  - cluster_number_scaling - The clustering algorithm tries to determine the optimal number of clusters (i.e. the number of peptide pairs in the LC-MS data set). If neighboring clusters appear in the same color, the cluster number is too low. If a single cluster contains two colors, the cluster number is too high. The cluster number can be adjusted by this scaling factor.
  - optimal_silhouette_tolerance - The clustering algorithm tries to maximize the average-silhouette-width, see details in reference. The parameter specifies the relative tolerance (in %) by which the optimum can deviate from the maximum.

  <b>References:</b>
  @n L. Nilse, M. Sturm, D. Trudgian, M. Salek, P. Sims, K. Carroll, S. Hubbard, "SILACAnalyzer - a tool for differential quantitation of stable isotope derived data", unpublished.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_SILACAnalyzer.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

struct SILACData
{
  DoubleReal rt; // retention time
  DoubleReal mz; // m/Z mass charge ratio
  DoubleReal int1; // intensity at RT and m/z
  DoubleReal int2; // intensity at RT and m/z + isotope_distance
  DoubleReal int3; // intensity at RT and m/z + 2*isotope_distance
  DoubleReal int4; // intensity at RT and m/z + envelope_distance_light_medium
  DoubleReal int5; // intensity at RT and m/z + envelope_distance_light_medium + isotope_distance
  DoubleReal int6; // intensity at RT and m/z + envelope_distance_light_medium + 2*isotope_distance
  DoubleReal int7; // intensity at RT and m/z + envelope_distance_light_heavy
  DoubleReal int8; // intensity at RT and m/z + envelope_distance_light_heavy + isotope_distance
  DoubleReal int9; // intensity at RT and m/z + envelope_distance_light_heavy + 2*isotope_distance
  Int cluster_id; // ID number of the cluster the data point belongs to
  Int cluster_size; // number of points in cluster 'cluster_id'

  SILACData();
  SILACData(const SILACData &);
  SILACData(DoubleReal rt_, DoubleReal mz_, DoubleReal int1_, DoubleReal int2_, DoubleReal int3_, DoubleReal int4_, DoubleReal int5_, DoubleReal int6_, DoubleReal int7_, DoubleReal int8_, DoubleReal int9_);
  ~SILACData(){};
  SILACData &operator=(const SILACData &rhs);
  int operator==(const SILACData &rhs) const;
  int operator<(const SILACData &rhs) const;
};

/// Default constructor
inline SILACData::SILACData()
{
  rt = 0;
  mz = 0;
  int1 = 0;
  int2 = 0;
  int3 = 0;
  int4 = 0;
  int5 = 0;
  int6 = 0;
  int7 = 0;
  int8 = 0;
  int9 = 0;
  cluster_id = 0;
  cluster_size = 0;
}
/// Copy constructor
inline SILACData::SILACData(const SILACData &copyin)
{
  rt = copyin.rt;
  mz = copyin.mz;
  int1 = copyin.int1;
  int2 = copyin.int2;
  int3 = copyin.int3;
  int4 = copyin.int4;
  int5 = copyin.int5;
  int6 = copyin.int6;
  int7 = copyin.int7;
  int8 = copyin.int8;
  int9 = copyin.int9;
  cluster_id = copyin.cluster_id;
  cluster_size = copyin.cluster_size;
}
/// Detailed constructor
inline SILACData::SILACData(DoubleReal rt_, DoubleReal mz_, DoubleReal int1_, DoubleReal int2_, DoubleReal int3_, DoubleReal int4_, DoubleReal int5_, DoubleReal int6_, DoubleReal int7_, DoubleReal int8_, DoubleReal int9_)
{
  rt = rt_;
  mz = mz_;
  int1 = int1_;
  int2 = int2_;
  int3 = int3_;
  int4 = int4_;
  int5 = int5_;
  int6 = int6_;
  int7 = int7_;
  int8 = int8_;
  int9 = int9_;
  cluster_id = 0;
  cluster_size = 0;
}

SILACData& SILACData::operator=(const SILACData &rhs)
{
  this->rt = rhs.rt;
  this->mz = rhs.mz;
  this->int1 = rhs.int1;
  this->int2 = rhs.int2;
  this->int3 = rhs.int3;
  this->int4 = rhs.int4;
  this->int5 = rhs.int5;
  this->int6 = rhs.int6;
  this->int7 = rhs.int7;
  this->int8 = rhs.int8;
  this->int9 = rhs.int9;
  this->cluster_id = rhs.cluster_id;
  this->cluster_size = rhs.cluster_size;
  return *this;
}

int SILACData::operator==(const SILACData &rhs) const
{
  if( this->rt != rhs.rt) return false;
  if( this->mz != rhs.mz) return false;
  if( this->int1 != rhs.int1) return false;
  if( this->int2 != rhs.int2) return false;
  if( this->int3 != rhs.int3) return false;
  if( this->int4 != rhs.int4) return false;
  if( this->int5 != rhs.int5) return false;
  if( this->int6 != rhs.int6) return false;
  if( this->int7 != rhs.int7) return false;
  if( this->int8 != rhs.int8) return false;
  if( this->int9 != rhs.int9) return false;
  if( this->cluster_id != rhs.cluster_id) return false;
  if( this->cluster_size != rhs.cluster_size) return false;
  return true;
}

int SILACData::operator<(const SILACData &rhs) const // required for built-in STL functions like sort
{
  if ( this->cluster_size == rhs.cluster_size && this->cluster_id < rhs.cluster_id) return true;
  if ( this->cluster_size < rhs.cluster_size ) return true;
  return false;
}

class TOPPSILACAnalyzer
      : public TOPPBase
{
  public:
    TOPPSILACAnalyzer()
        : TOPPBase("SILACAnalyzer","Determination of peak ratios in LC-MS data",true)
    {
    }

    void registerOptionsAndFlags_()
    {
	  registerStringOption_("type","<name>","","SILAC experiment type\n",true);
	  setValidStrings_("type", getToolList()[toolName_()] );

      registerInputFile_("in","<file>","","input file");
      setValidFormats_("in",StringList::create("mzML"));
      registerOutputFile_("out","<file>","","output file", false);
      setValidFormats_("out",StringList::create("consensusXML"));
      registerOutputFile_("out_visual","<file>","","output file containing cluster information",false);
      setValidFormats_("out_visual",StringList::create("featureXML"));

      registerFlag_("silac_debug","Enables writing of debug information",true);    
      registerSubsection_("algorithm","Algorithm parameters section");

    }

    Param getSubsectionDefaults_(const String& /*section*/) const
    {
      String type = getStringOption_("type");
      Int SILAC_type = (type=="double" ?  2 : 3 );
      Param tmp;

      if (SILAC_type == 2)
	  {
			tmp.setValue("mass_separation_light_heavy",6.0202, "mass gap between light and heavy isotopic envelopes, [Da]" );
	  } 
	  else if (SILAC_type == 3)
	  {
			tmp.setValue("mass_separation_light_medium" , 6.0202, "mass gap between light and medium isotopic envelopes, [Da]");
			tmp.setValue("mass_separation_light_heavy" , 8.0202, "mass gap between light and heavy isotopic envelopes, [Da]");
	  }
		
	  tmp.setValue("charge_min", 2, "Charge state range begin");
      tmp.setMinInt("charge_min", 1);

      tmp.setValue("charge_max", 3, "Charge state range end");
      tmp.setMinInt("charge_max",1);

      tmp.setValue("intensity_cutoff", 5000.0, "intensity cutoff");
      tmp.setMinFloat("intensity_cutoff", 0.0);

      tmp.setValue("mz_step_width", 0.01, "step width with which the (interpolated) spectrum is scanned, m/Z (Th)");
      tmp.setMinFloat("mz_step_width",0.0);

      tmp.setValue("rt_scaling",0.05,"scaling factor of retention times (Cluster height [s] an\ncluster width [Th] should be of the same order. The clustering algorithms work better for\nsymmetric clusters.)");
      tmp.setMinFloat("rt_scaling", 0.0);

      tmp.setValue("optimal_silhouette_tolerance",0.0,"The partition with most clusters is chosen, which deviates from the optimal silhouette width at most by this percentage.");
      tmp.setMinFloat("optimal_silhouette_tolerance",0.0);
      tmp.setMaxFloat("optimal_silhouette_tolerance",100.0);

      tmp.setValue("cluster_number_scaling", 1.0, "scaling factor of the number of clusters (The average-silhouette-width\nalgorithm returns an 'optimal' number of clusters. This number might need\nto be adjusted by this factor.)");
      tmp.setMinFloat("cluster_number_scaling", 0.0);

      tmp.setValue("cluster_min", 0, "Start of the clusters range to be plotted by the gnuplot script");
      tmp.setMinInt("cluster_min", 0);

      tmp.setValue("cluster_max", 2, "End of the clusters range to be plotted by the gnuplot script");
      tmp.setMinInt("cluster_max", 0);      


      return tmp;
    }

      ExitCodes main_(int , const char**)
    {
      //-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------
	  String type = getStringOption_("type");
		DoubleReal mass_separation_light_medium;
		if (type=="double") {
		  mass_separation_light_medium = getParam_().getValue("algorithm:mass_separation_light_heavy");
		}
		else {
		  mass_separation_light_medium = getParam_().getValue("algorithm:mass_separation_light_medium");
		}
	  DoubleReal mass_separation_light_heavy = getParam_().getValue("algorithm:mass_separation_light_heavy");

      UInt charge_min = getParam_().getValue("algorithm:charge_min");
      UInt charge_max = getParam_().getValue("algorithm:charge_max");

      DoubleReal mz_step_width = getParam_().getValue("algorithm:mz_step_width");
      DoubleReal intensity_cutoff = getParam_().getValue("algorithm:intensity_cutoff");
      DoubleReal rt_scaling = getParam_().getValue("algorithm:rt_scaling");
      DoubleReal optimal_silhouette_tolerance = getParam_().getValue("algorithm:optimal_silhouette_tolerance");
      DoubleReal cluster_number_scaling = getParam_().getValue("algorithm:cluster_number_scaling");
      int cluster_min = getParam_().getValue("algorithm:cluster_min");
      int cluster_max = getParam_().getValue("algorithm:cluster_max");

      String in = getStringOption_("in");
      String out = getStringOption_("out");
      String out_visual = getStringOption_("out_visual");

      //output variables
      ConsensusMap all_pairs;
	  all_pairs.getFileDescriptions()[0].filename = in;
	  all_pairs.getFileDescriptions()[0].label = "light";
	  all_pairs.getFileDescriptions()[1].filename = in;
	  all_pairs.getFileDescriptions()[1].label = "medium";
      all_pairs.getFileDescriptions()[2].filename = in;
      all_pairs.getFileDescriptions()[2].label = "heavy";
      all_pairs.setExperimentType("silac");
      FeatureMap<> all_cluster_points;

      //--------------------------------------------------------------
      // determine file name for debug output
      //--------------------------------------------------------------
      String debug_trunk = in;
      if (in.has('.'))
      {
        debug_trunk = in.substr(0,-SignedSize(in.suffix('.').length())-1);
      }

      // number of clusters found for each charge state (filled with best_n, need to remember for gnuplot script)
      int cluster_number[] = {1,1,1,1,1,1,1,1,1,1}; // maybe though with a vector even if more that 10 charge states are doubtfull?

      //iterate over all charge states
      for (UInt charge=charge_min; charge<=charge_max; ++charge)
      {
		std::cout << "charge state: " << charge << "+" << std::endl;
		DoubleReal isotope_distance = 1.0 / (DoubleReal)charge;
		DoubleReal envelope_distance_light_medium = mass_separation_light_medium / (DoubleReal)charge;
		DoubleReal envelope_distance_light_heavy = mass_separation_light_heavy / (DoubleReal)charge;
        // For each charge state run the experimental data (exp) are loaded again. Either the raw data (exp) or the distance matrix (distance_matrix) are in memory which keeps the memory footprint low.

        //-------------------------------------------------------------
        // loading input
        //-------------------------------------------------------------
        MzMLFile file;
        MSExperiment<Peak1D> exp;

        file.setLogType(log_type_);
        file.load(in,exp);

        //set input map size (only once)
        if (charge==charge_min)
        {
          exp.updateRanges();
		  all_pairs.getFileDescriptions()[0].size = exp.getSize();
		  all_pairs.getFileDescriptions()[1].size = exp.getSize();
          all_pairs.getFileDescriptions()[2].size = exp.getSize();
        }

        //-------------------------------------------------------------
        // build SILACData structure
        //-------------------------------------------------------------
        ProgressLogger logger_;
        std::vector<SILACData> data;

        logger_.setLogType(log_type_);
        logger_.startProgress(0,exp.size(),"reducing raw data");

        // scan over the entire experiment and write to data structure
        for (MSExperiment<>::Iterator rt_it=exp.begin(); rt_it!=exp.end(); ++rt_it)
        {
          logger_.setProgress(rt_it-exp.begin());
          Size number_data_points = rt_it->size();
          // spectra with less than 10 data points are being ignored
          if (number_data_points>=10) {
            // read one OpenMS spectrum into GSL structure
            std::vector<DoubleReal> mz_vec;
            std::vector<DoubleReal> intensity_vec;
            mz_vec.resize(number_data_points);
            intensity_vec.resize(number_data_points);
            Int j = 0;
            for (MSSpectrum<>::Iterator mz_it=rt_it->begin(); mz_it!=rt_it->end(); ++mz_it)
            {
              mz_vec[j] = mz_it->getMZ();
              intensity_vec[j] = mz_it->getIntensity();
              ++j;
            }
            DoubleReal mz_min = mz_vec[0];
            DoubleReal mz_max = mz_vec[number_data_points-1];
            // linear interpolation
            // used for the detection of pairs (spline overestimates at noise level)
            gsl_interp_accel *acc = gsl_interp_accel_alloc();
            gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, number_data_points);
            gsl_spline_init(spline, &*mz_vec.begin(), &*intensity_vec.begin(), number_data_points);
            // spline interpolation
            // used for exact ratio calculation (more accurate when real peak pairs are present)
            gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
            gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_cspline, number_data_points);
            gsl_spline_init(spline2, &*mz_vec.begin(), &*intensity_vec.begin(), number_data_points);
            for (DoubleReal mz=mz_min+isotope_distance; mz<mz_max-envelope_distance_light_heavy-3*isotope_distance; mz+=mz_step_width)
            {
              DoubleReal int_lin1 = gsl_spline_eval (spline, mz, acc);
              DoubleReal int_lin2 = gsl_spline_eval (spline, mz+isotope_distance, acc);
              DoubleReal int_lin3 = gsl_spline_eval (spline, mz+2*isotope_distance, acc);
			  DoubleReal int_lin4 = gsl_spline_eval (spline, mz+envelope_distance_light_medium, acc);
			  DoubleReal int_lin5 = gsl_spline_eval (spline, mz+envelope_distance_light_medium+isotope_distance, acc);
			  DoubleReal int_lin6 = gsl_spline_eval (spline, mz+envelope_distance_light_medium+2*isotope_distance, acc);
			  DoubleReal int_lin7 = gsl_spline_eval (spline, mz+envelope_distance_light_heavy, acc);
			  DoubleReal int_lin8 = gsl_spline_eval (spline, mz+envelope_distance_light_heavy+isotope_distance, acc);
			  DoubleReal int_lin9 = gsl_spline_eval (spline, mz+envelope_distance_light_heavy+2*isotope_distance, acc);
			  DoubleReal int_spline1 = gsl_spline_eval (spline2, mz, acc2);
              DoubleReal int_spline2 = gsl_spline_eval (spline2, mz+isotope_distance, acc2);
              DoubleReal int_spline3 = gsl_spline_eval (spline2, mz+2*isotope_distance, acc2);
			  DoubleReal int_spline4 = gsl_spline_eval (spline2, mz+envelope_distance_light_medium, acc2);
			  DoubleReal int_spline5 = gsl_spline_eval (spline2, mz+envelope_distance_light_medium+isotope_distance, acc2);
			  DoubleReal int_spline6 = gsl_spline_eval (spline2, mz+envelope_distance_light_medium+2*isotope_distance, acc2);
			  DoubleReal int_spline7 = gsl_spline_eval (spline2, mz+envelope_distance_light_heavy, acc2);
			  DoubleReal int_spline8 = gsl_spline_eval (spline2, mz+envelope_distance_light_heavy+isotope_distance, acc2);
			  DoubleReal int_spline9 = gsl_spline_eval (spline2, mz+envelope_distance_light_heavy+2*isotope_distance, acc2);
				
			  bool condDouble1 = (int_lin1 >= intensity_cutoff) && (int_lin2 >= intensity_cutoff) && (int_lin3 >= intensity_cutoff) && (int_lin7 >= intensity_cutoff) && (int_lin8 >= intensity_cutoff) && (int_lin9 >= intensity_cutoff); // all six intensities peak simultaneously
			  bool condDouble2 = (int_spline1 >= int_spline2) && (int_spline2 >= int_spline3) && (int_spline7 >= int_spline8) && (int_spline8 >= int_spline9); // isotopic peaks within one envelop decrease
			  bool condTriple1 = (int_lin1 >= intensity_cutoff) && (int_lin2 >= intensity_cutoff) && (int_lin3 >= intensity_cutoff) && (int_lin4 >= intensity_cutoff) && (int_lin5 >= intensity_cutoff) && (int_lin6 >= intensity_cutoff) && (int_lin7 >= intensity_cutoff) && (int_lin8 >= intensity_cutoff) && (int_lin9 >= intensity_cutoff); // all nine intensities peak simultaneously
			  bool condTriple2 = (int_spline1 >= int_spline2) && (int_spline2 >= int_spline3) && (int_spline4 >= int_spline5) && (int_spline5 >= int_spline6) && (int_spline7 >= int_spline8) && (int_spline8 >= int_spline9); // isotopic peaks within one envelop decrease
              if ((type=="double" && condDouble1 && condDouble2) || (type=="triple" && condTriple1 && condTriple2))
              {
                data.push_back(SILACData(rt_it->getRT(),mz,int_spline1,int_spline2,int_spline3,int_spline4,int_spline5,int_spline6,int_spline7,int_spline8,int_spline9));
              }
            }

            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            gsl_spline_free(spline2);
            gsl_interp_accel_free(acc2);
          }
        }
        exp.clear(true);
        logger_.endProgress();

        //-------------------------------------------------------------
        // generate distance matrix and copy distance matrix
        //-------------------------------------------------------------
        DistanceMatrix<Real> distance_matrix;
        distance_matrix.resize(data.size(),1);
        for (Size i=0; i < data.size(); ++i)
        {
          for (Size j=0; j < i; ++j)
          {
            // shrink RT by factor rt_scaling in order to make clusters more symmetric
            distance_matrix.setValueQuick(i,j, sqrt( (data[i].rt-data[j].rt)*(data[i].rt-data[j].rt)*rt_scaling*rt_scaling+(data[i].mz-data[j].mz)*(data[i].mz-data[j].mz) ) );
          }
        }

        DistanceMatrix<Real> distance_matrix_copy(distance_matrix); //clustering method will mess with input matrix

        //-------------------------------------------------------------
        // conduct clustering
        //-------------------------------------------------------------
        //~ ClusterHierarchical ch;
		AverageLinkage al;
        al.setLogType(log_type_);
        std::vector< BinaryTreeNode > tree;
        al(distance_matrix_copy, tree, std::numeric_limits<float>::max());

        //-----------------------------------------------------------------
        // find number of clusters which maximizes average silhouette width
        //-----------------------------------------------------------------

        ClusterAnalyzer ca;
        //choose asw that deviates at most the given percentage (optimal_silhouette_tolerance) from the max asw and contains the most clusters
        std::vector< Real >asw = ca.averageSilhouetteWidth(tree,distance_matrix);
        std::vector< Real >::iterator max_el(max_element(asw.begin(),asw.end()));
        //~ std::vector< Real >::iterator max_el(max_element((asw.end()-((Int)data.size()/10) ),asw.end()));//only the first size/10 steps are reviewed
        int best_n = (int)tree.size();
        Real max_deviation((*max_el)*(optimal_silhouette_tolerance/100));
        for (Size i = 0; i < asw.size(); ++i)
        {
          if(std::fabs(asw[i]-(*max_el))<=max_deviation)
          {
            best_n = int(tree.size() - i);
            break;
          }
        }

        //-------------------------------------------------------------
        // choose appropriate(best) partition of data from best_n
        //-------------------------------------------------------------
        best_n = int(cluster_number_scaling * best_n); // slightly increase cluster number
        std::vector< std::vector<Size> > best_n_clusters;
        ca.cut(best_n,best_n_clusters,tree);
        cluster_number[charge] = best_n;

        //-------------------------------------------------------------
        // count data points in each cluster
        //-------------------------------------------------------------
        std::vector<SignedSize> cluster_size(best_n,0); // number of points per cluster
        for (Size i=0; i < cluster_size.size(); ++i)
        {
          cluster_size[i] = best_n_clusters[i].size();
        }

        //--------------------------------------------------------------
        // fill in cluster_id and cluster_size in SILACData structure
        //--------------------------------------------------------------
        for (Size i=0; i < best_n_clusters.size(); ++i)
        {
          for (Size j=0; j < best_n_clusters[i].size(); ++j)
          {
            data[best_n_clusters[i][j]].cluster_id = (Int)i;
            data[best_n_clusters[i][j]].cluster_size =(Int)best_n_clusters[i].size();
          }
        }
        std::sort(data.begin(),data.end());
        std::reverse(data.begin(),data.end()); // largest clusters first

        //--------------------------------------------------------------
        // update cluster_id
        //--------------------------------------------------------------
        Int k = -1;
        Int new_id = best_n-1;
        for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
        {
          if (it->cluster_id != k) ++new_id;
          k = it->cluster_id;
          it->cluster_id = new_id;
        }
        for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
        {
          it->cluster_id = it->cluster_id - best_n;
        }

        //--------------------------------------------------------------
        // update cluster_size
        //--------------------------------------------------------------
        cluster_size = std::vector<SignedSize>(best_n,0);
        k = -1;
        for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
        {
          ++cluster_size[it->cluster_id];
        }

        //--------------------------------------------------------------
        //create consensus features
        //--------------------------------------------------------------
        if (out!="")
        {
          for (int i=0; i<best_n; ++i)
          {
            DoubleReal rt = 0.0;
            DoubleReal mz = 0.0;
			DoubleReal int_l = 0.0;
			DoubleReal int_m = 0.0;
            DoubleReal int_h = 0.0;
            // intensity vectors used for linear regression
			std::vector<DoubleReal> i1(3*cluster_size[i]);
			std::vector<DoubleReal> i2(3*cluster_size[i]);
            std::vector<DoubleReal> i3(3*cluster_size[i]);
            UInt j=0;
            for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
            {
              if (it->cluster_id==i)
              {
				i1[3*j] = it->int1;
				i2[3*j] = it->int4;
                i3[3*j] = it->int7;
				i1[3*j+1] = it->int2;
				i2[3*j+1] = it->int5;
                i3[3*j+1] = it->int8;
				i1[3*j+2] = it->int3;
				i2[3*j+2] = it->int6;
                i3[3*j+2] = it->int9;

                rt += it->rt;
                if (it->int1>int_l)
                {
                  int_l = it->int1;
                  mz = it->mz; // m/z of maximum intensity in the light cluster
                }
                if (it->int2>int_l)
                {
                  int_l = it->int2;
                  mz = it->mz + isotope_distance;
                }
                if (it->int3>int_l)
                {
                  int_l = it->int3;
                  mz = it->mz + 2.0 * isotope_distance;
                }
				if (it->int4>int_m)
				{
					int_m = it->int4;
				}
				if (it->int5>int_m)
				{
					int_m = it->int5;
				}
				if (it->int6>int_m)
				{
					int_m = it->int6;
				}
				if (it->int7>int_h)
				{
					int_h = it->int7;
				}
				if (it->int8>int_h)
				{
					int_h = it->int8;
				}
				if (it->int9>int_h)
				{
					int_h = it->int9;
				}
				++j;
              }
            }
            rt /= (DoubleReal)(cluster_size[i]); // average retention time
			Math::LinearRegression linear_reg_light_medium;
			Math::LinearRegression linear_reg_light_heavy;
			linear_reg_light_medium.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i2.begin());
			linear_reg_light_heavy.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i3.begin());
            //create consensus feature for light medium SILAC pair
            ConsensusFeature pair_light_medium;
            pair_light_medium.setRT(rt);
            pair_light_medium.setMZ(mz);
            pair_light_medium.setIntensity(linear_reg_light_medium.getSlope());
            pair_light_medium.setCharge(charge);
            pair_light_medium.setQuality(linear_reg_light_medium.getRSquared());
            FeatureHandle handle;
            handle.setRT(rt);
            handle.setMZ(mz);
            handle.setIntensity(int_l);
            handle.setCharge(charge);
            handle.setMapIndex(0);
            handle.setUniqueId(i);
            pair_light_medium.insert(handle);
            handle.setRT(rt);
            handle.setMZ(mz+envelope_distance_light_medium);
            handle.setIntensity(int_m);
            handle.setCharge(charge);
            handle.setMapIndex(1);
            handle.setUniqueId(i);
            pair_light_medium.insert(handle);
			if (type=="triple") {
			  all_pairs.push_back(pair_light_medium);
			}
			//create consensus feature for light heavy SILAC pair
            ConsensusFeature pair_light_heavy;
			pair_light_heavy.setRT(rt);
			pair_light_heavy.setMZ(mz);
			pair_light_heavy.setIntensity(linear_reg_light_heavy.getSlope());
			pair_light_heavy.setCharge(charge);
			pair_light_heavy.setQuality(linear_reg_light_heavy.getRSquared());
			handle.setRT(rt);
			handle.setMZ(mz);
			handle.setIntensity(int_l);
			handle.setCharge(charge);
			handle.setMapIndex(0);
			handle.setUniqueId(i);
			pair_light_heavy.insert(handle);
			handle.setRT(rt);
			handle.setMZ(mz+envelope_distance_light_heavy);
			handle.setIntensity(int_h);
			handle.setCharge(charge);
			handle.setMapIndex(2);
			handle.setUniqueId(i);
			pair_light_heavy.insert(handle);
			all_pairs.push_back(pair_light_heavy);
          }
        }


        //--------------------------------------------------------------
        //create features (for visualization)
        //--------------------------------------------------------------
        if (out_visual!="")
        {
          std::vector<String> colors;
          // 16 HTML colors
          colors.push_back("#00FFFF");
          colors.push_back("#000000");
          colors.push_back("#0000FF");
          colors.push_back("#FF00FF");
          colors.push_back("#008000");
          colors.push_back("#808080");
          colors.push_back("#00FF00");
          colors.push_back("#800000");
          colors.push_back("#000080");
          colors.push_back("#808000");
          colors.push_back("#800080");
          colors.push_back("#FF0000");
          colors.push_back("#C0C0C0");
          colors.push_back("#008080");
          colors.push_back("#FFFF00");

		  for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
          {
            //visualize the light variant
            Feature cluster_point;
            cluster_point.setRT(it->rt);
            cluster_point.setMZ(it->mz);
            DoubleReal intensity = std::max(it->int1,it->int2);
            intensity = std::max(intensity,it->int3);
            cluster_point.setIntensity(intensity);
            cluster_point.setCharge(charge);
            cluster_point.setMetaValue("cluster_id",it->cluster_id);
            cluster_point.setMetaValue("color",colors[it->cluster_id%colors.size()]);
            all_cluster_points.push_back(cluster_point);
          }
          // required, as somehow the order of features on some datasets between Win & Linux is different and thus the TOPPtest might fail
          all_cluster_points.sortByPosition();
        }





		//-------------------------------------------------------------
		// generate debug output
		//-------------------------------------------------------------
		// strings repeatedly used in debug output
		String light_medium_string = String(0.01*floor(mass_separation_light_medium*100+0.5)); 
		String light_heavy_string = String(0.01*floor(mass_separation_light_heavy*100+0.5));
		  
		if (getFlag_("silac_debug"))
        {
		  String debug_suffix;
		  if (type=="double") {
			debug_suffix = "_" + light_heavy_string + "Da_" + String(charge) +"+";
		  }
		  else {
			debug_suffix = "_" + light_medium_string + "Da_" + light_heavy_string + "Da_" + String(charge) +"+";
		  }
          // names of dat files
          String debug_dat = debug_trunk + debug_suffix + ".dat";
          String debug_clusters_dat = debug_trunk + debug_suffix + "_clusters.dat";

          // write all cluster data points to *_clusters.dat
          std::ofstream stream_clusters(debug_clusters_dat.c_str());
		  if (type=="double") {
			stream_clusters << "cluster_id cluster_size rt mz int1 int2 int3 int7 int8 int9" << std::endl;
		  }
		  else {
			stream_clusters << "cluster_id cluster_size rt mz int1 int2 int3 int4 int5 int6 int7 int8 int9" << std::endl;
		  }
          Int current_id = -1;
          for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
          {
            if (it->cluster_id != current_id) stream_clusters << std::endl << std::endl;
			if (type=="double") {
			  stream_clusters << it->cluster_id << " " << it->cluster_size << " "<< it->rt << " " << it->mz << " " << it->int1 << " " << it->int2 << " " << it->int3 << " " << it->int7 << " " << it->int8 << " " << it->int9 << std::endl;
			}
			else {
			  stream_clusters << it->cluster_id << " " << it->cluster_size << " "<< it->rt << " " << it->mz << " " << it->int1 << " " << it->int2 << " " << it->int3 << " " << it->int4 << " " << it->int5 << " " << it->int6 << " " << it->int7 << " " << it->int8 << " " << it->int9 << std::endl;
			}
            current_id = it->cluster_id;
          }
          stream_clusters.close();

          // write ratios of all cluster to *.dat
          std::ofstream stream_ratios(debug_dat.c_str());
		  if (type=="double") {
			stream_ratios << "cluster_id cluster_size rt mz ratio_light_heavy" << std::endl;
		  }
		  else {
			stream_ratios << "cluster_id cluster_size rt mz ratio_light_medium ratio_light_heavy" << std::endl;
		  }
          for (int i=0; i<best_n;++i)
          {
            DoubleReal rt = 0.0;
            DoubleReal mz = 0.0;
			DoubleReal int_l = 0.0;
			DoubleReal int_m = 0.0;
            DoubleReal int_h = 0.0;
            // intensity vectors used for linear regression
			std::vector<DoubleReal> i1(3*cluster_size[i]);
			std::vector<DoubleReal> i2(3*cluster_size[i]);
            std::vector<DoubleReal> i3(3*cluster_size[i]);
            UInt j=0;
             for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
            {
              if (it->cluster_id==i)
              {
				i1[3*j] = it->int1;
				i2[3*j] = it->int4;
                i3[3*j] = it->int7;
				i1[3*j+1] = it->int2;
				i2[3*j+1] = it->int5;
                i3[3*j+1] = it->int8;
				i1[3*j+2] = it->int3;
				i2[3*j+2] = it->int6;
                i3[3*j+2] = it->int9;

                rt += it->rt;
                if (it->int1>int_l)
                {
                  int_l = it->int1;
                  mz = it->mz;
                }
                if (it->int2>int_l)
                {
                  int_l = it->int2;
                  mz = it->mz + isotope_distance;
                }
                if (it->int3>int_l)
                {
                  int_l = it->int3;
                  mz = it->mz + 2.0 * isotope_distance;
                }
				if (it->int4>int_m)
				{
				  int_m = it->int4;
				}
				if (it->int5>int_m)
				{
				  int_m = it->int5;
				}
				if (it->int6>int_m)
				{
				  int_m = it->int6;
				}
				if (it->int7>int_h)
				{
				  int_h = it->int7;
				}
				if (it->int8>int_h)
				{
				  int_h = it->int8;
				}
				if (it->int9>int_h)
				{
				  int_h = it->int9;
				}
				++j;
              }
            }
            rt /= (DoubleReal)(cluster_size[i]); // average retention time
			Math::LinearRegression linear_reg_light_medium;
			Math::LinearRegression linear_reg_light_heavy;
			linear_reg_light_medium.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i2.begin());
			linear_reg_light_heavy.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i3.begin());
			if (type=="double") {
			  stream_ratios << i << " " << cluster_size[i] << " " << rt << " " << mz << " " << linear_reg_light_heavy.getSlope() << std::endl;
			}
			else {
			  stream_ratios << i << " " << cluster_size[i] << " " << rt << " " << mz << " " << linear_reg_light_medium.getSlope() << " " << linear_reg_light_heavy.getSlope() << std::endl;
			}
          }
          stream_ratios.close();
        }

      } //end iterate over all charge states


      //--------------------------------------------------------------
      //Store output
      //--------------------------------------------------------------
      if (out!="")
      {

        // assign unique ids
        all_pairs.applyMemberFunction(&UniqueIdInterface::setUniqueId);

        //annotate output with data processing info
        addDataProcessing_(all_pairs, getProcessingInfo_(DataProcessing::QUANTITATION));

        ConsensusXMLFile c_file;
        c_file.store(out,all_pairs);
      }

      if (out_visual!="")
      {
        // assign unique ids
        all_cluster_points.applyMemberFunction(&UniqueIdInterface::setUniqueId);

        FeatureXMLFile f_file;
        f_file.store(out_visual,all_cluster_points);
      }

      //--------------------------------------------------------------
      //write gnuplot script
      //--------------------------------------------------------------
	  // strings repeatedly used in debug output
	  String light_medium_string = String(0.01*floor(mass_separation_light_medium*100+0.5)); 
	  String light_heavy_string = String(0.01*floor(mass_separation_light_heavy*100+0.5));
	  String rt_scaling_string = String(0.01*floor(rt_scaling*100+0.5));
	  String optimal_silhouette_tolerance_string = String(0.01*floor(optimal_silhouette_tolerance*100+0.5));
	  String cluster_number_scaling_string = String(0.01*floor(cluster_number_scaling*100+0.5));
		
	  if (getFlag_("silac_debug"))
      {
        // first lines of the gnuplot script
        String debug_gnuplotscript = debug_trunk + ".input";
        std::ofstream stream_gnuplotscript(debug_gnuplotscript.c_str());
        stream_gnuplotscript << "set terminal postscript eps enhanced colour" << std::endl;
        stream_gnuplotscript << "set size 2.0, 2.0" << std::endl;
        stream_gnuplotscript << "set size square" << std::endl << std::endl;
        //iterate over all charge states
        for (Size charge=charge_min; charge<=charge_max; ++charge)
        {
		  String debug_suffix;
		  if (type=="double") {
			debug_suffix = "_" + light_heavy_string + "Da_" + String(charge) +"+";
		  }
		  else {
			debug_suffix = "_" + light_medium_string + "Da_" + light_heavy_string + "Da_" + String(charge) +"+";
		  }
          // names of dat files
          String debug_dat = debug_trunk + debug_suffix + ".dat";
          String debug_clusters_dat = debug_trunk + debug_suffix + "_clusters.dat";
          // names of postscript files
		  String debug_ratios_light_medium = debug_trunk + debug_suffix + "_ratios_light_medium.eps";
		  String debug_ratios_light_heavy = debug_trunk + debug_suffix + "_ratios_light_heavy.eps";
          String debug_sizes = debug_trunk + debug_suffix + "_sizes.eps";
          String debug_clusters = debug_trunk + debug_suffix + "_clusters.eps";
		  String debug_clustersIntLightMedium = debug_trunk + debug_suffix + "_clustersIntLightMedium.eps";
		  String debug_clustersIntLightHeavy = debug_trunk + debug_suffix + "_clustersIntLightHeavy.eps";

          // write *_clusters.eps
          stream_gnuplotscript << "set output \"" + debug_clusters + "\"" << std::endl;
		  if (type=="double") {
			stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << "\\n mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
		  }
		  else {
			stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
		  }
          stream_gnuplotscript << "set xlabel \'m/Z (Th)\'" << std::endl;
          stream_gnuplotscript << "set ylabel \'RT (s)\'" << std::endl;
          stream_gnuplotscript << "plot";
          for (int i = cluster_min; i <= cluster_max; i++)// iterate over clusters
          {
            if (i != cluster_min) stream_gnuplotscript << ",";
            stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 4:3 title \"cluster " + String(i) + "\"";
          }
          stream_gnuplotscript << std::endl;

		  // write *_clustersIntLightMedium.eps
		  if (type=="triple") {
		    stream_gnuplotscript << "set output \"" + debug_clustersIntLightMedium + "\"" << std::endl;
		    stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
		    stream_gnuplotscript << "set xlabel \'intensity at m/Z\'" << std::endl;
		    stream_gnuplotscript << "set ylabel \'intensity at m/Z + " + light_medium_string + "Th\'" << std::endl;
		    stream_gnuplotscript << "plot";
		    for (int i = cluster_min; i <= cluster_max; i++)// iterate over clusters
		    {
			  if (i != cluster_min) stream_gnuplotscript << ",";
			  stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 5:8 title \"cluster " + String(i) + "\"";
		    }
		    stream_gnuplotscript << std::endl;
		  }
			
		  // write *_clustersIntLightHeavy.eps
		  stream_gnuplotscript << "set output \"" + debug_clustersIntLightHeavy + "\"" << std::endl;
		  if (type=="double") {
			stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << "\\n mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
		  }
		  else {
			stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
		  }
		  stream_gnuplotscript << "set xlabel \'intensity at m/Z\'" << std::endl;
		  stream_gnuplotscript << "set ylabel \'intensity at m/Z + " + light_heavy_string + "Th\'" << std::endl;
		  stream_gnuplotscript << "plot";
		  for (int i = cluster_min; i <= cluster_max; i++)// iterate over clusters
		  {
			if (i != cluster_min) stream_gnuplotscript << ",";
			if (type=="double") {
				stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 5:8 title \"cluster " + String(i) + "\"";
			}
			else {
				stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 5:11 title \"cluster " + String(i) + "\"";
			}
		  }
		  stream_gnuplotscript << std::endl;
			
		  // write *_ratios_light_medium.eps
		  if (type=="triple") {
		    stream_gnuplotscript << "set output \"" + debug_ratios_light_medium + "\"" << std::endl;
			stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
		    stream_gnuplotscript << "set nokey" << std::endl;
		    stream_gnuplotscript << "set xlabel \'m/Z\'" << std::endl;
		    stream_gnuplotscript << "set ylabel \'ratio light medium\'" << std::endl;
			stream_gnuplotscript << "plot \'" + debug_dat + "\' using 4:5";
		    stream_gnuplotscript << std::endl;
		  }
			
		  // write *_ratios_light_heavy.eps
          stream_gnuplotscript << "set output \"" + debug_ratios_light_heavy + "\"" << std::endl;
		  if (type=="double") {
			stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << "\\n mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
		  }
		  else {
			stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
		  }
		  stream_gnuplotscript << "set nokey" << std::endl;
          stream_gnuplotscript << "set xlabel \'m/Z\'" << std::endl;
          stream_gnuplotscript << "set ylabel \'ratio light heavy\'" << std::endl;
		  if (type=="double") {
			stream_gnuplotscript << "plot \'" + debug_dat + "\' using 4:5";
		  }
		  else {
			stream_gnuplotscript << "plot \'" + debug_dat + "\' using 4:6";
		  }
          stream_gnuplotscript << std::endl;

          // write *_sizes.eps
          stream_gnuplotscript << "set output \"" + debug_sizes + "\"" << std::endl;
		  if (type=="double") {
			stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << "\\n mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
		  }
		  else {
			stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << "\\n mass separation light medium= " << light_medium_string << ", mass separation light heavy= " << light_heavy_string << " Da, charge = " << charge << "+\\n intensity cutoff = " << intensity_cutoff << ", rt scaling = " << rt_scaling_string << ", optimal silhouette tolerance = " << optimal_silhouette_tolerance_string << ", cluster number scaling = " << cluster_number_scaling_string << "\"" << std::endl;
		  }
		  stream_gnuplotscript << "set nokey" << std::endl;
          stream_gnuplotscript << "set xlabel \'cluster ID\'" << std::endl;
          stream_gnuplotscript << "set ylabel \'cluster size\'" << std::endl;
          stream_gnuplotscript << "plot \'" + debug_dat + "\' using 1:2";
          stream_gnuplotscript << std::endl;
        }
        stream_gnuplotscript.close();
      }


      return EXECUTION_OK;
    }
};

//@endcond

int main(int argc, const char** argv )
{
  TOPPSILACAnalyzer tool;
  return tool.main(argc,argv);
}
