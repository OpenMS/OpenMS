// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

//OpenMS includes
#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

//clustering
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/COMPARISON/CLUSTERING/AverageLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
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

	@brief Determines the ratio of peak pairs in LC-MS data.

	(1) data reduction
	(2) hierarchical clustering in RT-m/Z plane, determine cluster number by maximizing the average silhouette width
	(3) determine intensity ratios by linear regression for each cluster

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
	DoubleReal int4; // intensity at RT and m/z + envelope_distance
	DoubleReal int5; // intensity at RT and m/z + envelope_distance + isotope_distance
	DoubleReal int6; // intensity at RT and m/z + envelope_distance + 2*isotope_distance
	Int cluster_id; // ID number of the cluster the data point belongs to
	Int cluster_size; // number of points in cluster 'cluster_id'

	SILACData();
	SILACData(const SILACData &);
	SILACData(DoubleReal rt_, DoubleReal mz_, DoubleReal int1_, DoubleReal int2_, DoubleReal int3_, DoubleReal int4_, DoubleReal int5_, DoubleReal int6_);
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
	cluster_id = copyin.cluster_id;
	cluster_size = copyin.cluster_size;
}
/// Detailed constructor
inline SILACData::SILACData(DoubleReal rt_, DoubleReal mz_, DoubleReal int1_, DoubleReal int2_, DoubleReal int3_, DoubleReal int4_, DoubleReal int5_, DoubleReal int6_)
{
	rt = rt_;
	mz = mz_;
	int1 = int1_;
	int2 = int2_;
	int3 = int3_;
	int4 = int4_;
	int5 = int5_;
	int6 = int6_;
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
				: TOPPBase("SILACAnalyzer","Determination of peak ratios in LC-MS data",true,"0.6.2")
		{
		}

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","input file");
			setValidFormats_("in",StringList::create("mzData"));
			registerOutputFile_("out","<file>","","output file", false);
			setValidFormats_("out",StringList::create("consensusXML"));
			registerOutputFile_("out_visual","<file>","","output file containing cluster information",false);
			setValidFormats_("out_visual",StringList::create("featureXML"));

			registerFlag_("silac_debug","Enables writing of debug information",true);

			registerDoubleOption_("mass_separation","<dist>",6.0202,"m/z gap between light and heavy isotopic envelopes, [Da]",false);
			registerIntOption_("charge_min","<min>",2,"Charge state range begin",false);
			setMinInt_("charge_min",1);
			registerIntOption_("charge_max","<max>",3,"Charge state range end",false);
			setMinInt_("charge_max",1);
			registerDoubleOption_("intensity_cutoff","<double>",5000,"intensity cutoff",false,true);
			setMinFloat_("intensity_cutoff",0.0);
			registerDoubleOption_("mz_step_width","<double>",0.01,"step width with which the (interpolated) spectrum is scanned, m/Z (Th)",false,true);
			setMinFloat_("mz_step_width",0.0);
			registerDoubleOption_("rt_scaling","<double>",0.05,"scaling factor of retention times (Cluster height [s] an\ncluster width [Th] should be of the same order. The clustering algorithms work better for\nsymmetric clusters.)",false,true);
			setMinFloat_("rt_scaling",0.0);
			registerDoubleOption_("optimal_silhouette_tolerance","<double>",0.0,"The partition with most clusters is chosen, which deviates from the optimal silhouette width at most by this percentage.",false,true);
			setMinFloat_("optimal_silhouette_tolerance",0.0);
			setMaxFloat_("optimal_silhouette_tolerance",100.0);
			registerDoubleOption_("cluster_number_scaling","<double>",1.0,"scaling factor of the number of clusters (The average-silhouette-width\nalgorithm returns an 'optimal' number of clusters. This number might need\nto be adjusted by this factor.)",false,true); //still needed with blurred partition chooser?
			setMinFloat_("cluster_number_scaling",0.0);
			registerIntOption_("cluster_min","<min>",0,"Start of the clusters range to be plotted by the gnuplot script", false,true);
			setMinInt_("cluster_min",0);
			registerIntOption_("cluster_max","<max>",2,"End of the clusters range to be plotted by the gnuplot script", false,true);
			setMinInt_("cluster_max",0);
		}

			ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
			DoubleReal mass_separation = getDoubleOption_("mass_separation");

			UInt charge_min = getIntOption_("charge_min");
			UInt charge_max = getIntOption_("charge_max");

			DoubleReal mz_step_width = getDoubleOption_("mz_step_width");
			DoubleReal intensity_cutoff = getDoubleOption_("intensity_cutoff");
			DoubleReal rt_scaling = getDoubleOption_("rt_scaling");
			DoubleReal optimal_silhouette_tolerance = getDoubleOption_("optimal_silhouette_tolerance");
			DoubleReal cluster_number_scaling = getDoubleOption_("cluster_number_scaling");
			int cluster_min = getIntOption_("cluster_min");
			int cluster_max = getIntOption_("cluster_max");

			String in = getStringOption_("in");
			String out = getStringOption_("out");
			String out_visual = getStringOption_("out_visual");

			//output variables
			ConsensusMap all_pairs;
			all_pairs.getFileDescriptions()[0].filename = in;
			all_pairs.getFileDescriptions()[0].label = "light";
			all_pairs.getFileDescriptions()[0].size = 0;
			all_pairs.getFileDescriptions()[1].filename = in;
			all_pairs.getFileDescriptions()[1].label = "heavy";
			all_pairs.getFileDescriptions()[1].size = 0;
			all_pairs.setExperimentType("silac");
			FeatureMap<> all_cluster_points;

			//--------------------------------------------------------------
			// determine file name for debug output
			//--------------------------------------------------------------
			String debug_trunk = in;
			//std::cout << "in=" << in << std::endl;
			if (in.has('.'))
			{
				//std::cout << "yes" << std::endl;
				debug_trunk = in.substr(0,-in.suffix('.').length()-1);
				//std::cout << debug_trunk << std::endl;
			}

			// number of clusters found for each charge state (filled with best_n, need to remember for gnuplot script)
			int cluster_number[] = {1,1,1,1,1,1,1,1,1,1}; // maybe though with a vector even if more that 10 charge states are doubtfull?

			//iterate over all charge states
			for (UInt charge=charge_min; charge<=charge_max; ++charge)
			{
				std::cout << std::endl << "charge state: " << charge << "+" << std::endl;
				DoubleReal isotope_distance = 1.0 / (DoubleReal)charge;
				DoubleReal envelope_distance = mass_separation / (DoubleReal)charge;
				// For each charge state run the experimental data (exp) are loaded again. Either the raw data (exp) or the distance matrix (distance_matrix) are in memory which keeps the memory footprint low.

				//-------------------------------------------------------------
				// loading input
				//-------------------------------------------------------------
				MzDataFile file;
				MSExperiment<Peak1D> exp;

				file.setLogType(log_type_);
				file.load(in,exp);

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
					Int number_data_points = rt_it->size();
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
					  for (DoubleReal mz=mz_min+isotope_distance; mz<mz_max-envelope_distance-3*isotope_distance; mz+=mz_step_width)
					  {
						  DoubleReal int_lin1 = gsl_spline_eval (spline, mz, acc);
						  DoubleReal int_lin2 = gsl_spline_eval (spline, mz+envelope_distance, acc);
						  DoubleReal int_lin3 = gsl_spline_eval (spline, mz+isotope_distance, acc);
  						DoubleReal int_lin4 = gsl_spline_eval (spline, mz+envelope_distance+isotope_distance, acc);
	  					DoubleReal int_lin5 = gsl_spline_eval (spline, mz+2*isotope_distance, acc);
		  				DoubleReal int_lin6 = gsl_spline_eval (spline, mz+envelope_distance+2*isotope_distance, acc);
			  			DoubleReal int_spline1 = gsl_spline_eval (spline2, mz, acc2);
				  		DoubleReal int_spline2 = gsl_spline_eval (spline2, mz+envelope_distance, acc2);
					  	DoubleReal int_spline3 = gsl_spline_eval (spline2, mz+isotope_distance, acc2);
						  DoubleReal int_spline4 = gsl_spline_eval (spline2, mz+envelope_distance+isotope_distance, acc2);
  						DoubleReal int_spline5 = gsl_spline_eval (spline2, mz+2*isotope_distance, acc2);
	  					DoubleReal int_spline6 = gsl_spline_eval (spline2, mz+envelope_distance+2*isotope_distance, acc2);

		  				bool cond1 = (int_lin1 >= intensity_cutoff) && (int_lin2 >= intensity_cutoff) && (int_lin3 >= intensity_cutoff) && (int_lin4 >= intensity_cutoff) && (int_lin5 >= intensity_cutoff) && (int_lin6 >= intensity_cutoff); // all six intensities peak simultaneously
			  			bool cond2 = (int_spline3 <= int_spline1) && (int_spline5 <= int_spline3) && (int_spline4 <= int_spline2) && (int_spline6 <= int_spline4); // isotopic peaks within one envelop decrease
				  		if (cond1 && cond2)
					  	{
			  		  	data.push_back(SILACData(rt_it->getRT(),mz,int_spline1,int_spline3,int_spline5,int_spline2,int_spline4,int_spline6));
						  }
					  }

					  gsl_spline_free(spline);
  					gsl_interp_accel_free(acc);
	  				gsl_spline_free(spline2);
		  			gsl_interp_accel_free(acc2);
					}
				}
				exp.clear();
				logger_.endProgress();

				//-------------------------------------------------------------
				// generate distance matrix and copy distance matrix
				//-------------------------------------------------------------
				DistanceMatrix<Real> distance_matrix;
				distance_matrix.resize(data.size(),1);
				for(UInt i=0; i < data.size(); ++i)
				{
					for(UInt j=0; j < i; ++j)
					{
						// shrink RT by factor rt_scaling in order to make clusters more symmetric
						distance_matrix.setValueQuick(i,j, sqrt( (data[i].rt-data[j].rt)*(data[i].rt-data[j].rt)*rt_scaling*rt_scaling+(data[i].mz-data[j].mz)*(data[i].mz-data[j].mz) ) );
					}
				}

				DistanceMatrix<Real> distance_matrix_copy(distance_matrix); //clustering method will mess with input matrix

				//-------------------------------------------------------------
				// conduct clustering
				//-------------------------------------------------------------
				ClusterHierarchical ch;
				AverageLinkage al;
				std::vector< BinaryTreeNode > tree;
				ClusterAnalyzer ca;
				al.cluster(distance_matrix_copy,tree,DBL_MAX);

				//-----------------------------------------------------------------
				// find number of clusters which maximizes average silhouette width
				//-----------------------------------------------------------------

				//choose asw that deviates at most the given percentage (optimal_silhouette_tolerance) from the max asw and contains the most clusters
				std::vector< Real >asw = ca.averageSilhouetteWidth(tree,distance_matrix);
				std::vector< Real >::iterator max_el(max_element(asw.begin(),asw.end()));
				//~ std::vector< Real >::iterator max_el(max_element((asw.end()-((Int)data.size()/10) ),asw.end()));//only the first size/10 steps are reviewed
				size_t best_n = tree.size();
				Real max_deviation((*max_el)*(optimal_silhouette_tolerance/100));
				for(UInt i = 0; i < asw.size(); ++i)
				{
					if(std::fabs(asw[i]-(*max_el))<=max_deviation)
					{
						best_n = tree.size() - i;
						break;
					}
				}

				//~ debug output
				//~ for(UInt j = 0; j < asw.size(); ++j)
				//~ {
					//~ std::cout << "asw(" << j << "): " << asw[j] << std::endl;
				//~ }
				//~ std::cout << "best_n: " << best_n << std::endl;

				//~ dunn indices nogood here, if not only the the last half considered
				//~ std::vector<Real>dunn = ca.dunnIndices(tree,distance_matrix);
				//~ for(UInt j = 0; j < dunn.size(); ++j)
				//~ {
					//~ std::cout << "dunn: " << dunn[j] << std::endl;
				//~ }
				//~ max_el = (max_element(dunn.begin(),dunn.end()));
				//~ best_n = tree.size() - (max_el-dunn.begin());

				//-------------------------------------------------------------
				// choose appropriate(best) partition of data from best_n
				//-------------------------------------------------------------
				best_n = UInt(cluster_number_scaling * best_n); // slightly increase cluster number
				std::vector< std::vector<UInt> > best_n_clusters;
				ca.cut(best_n,best_n_clusters,tree);
				cluster_number[charge] = best_n;

				//-------------------------------------------------------------
				// count data points in each cluster
				//-------------------------------------------------------------
				std::vector<Int> cluster_size(best_n,0); // number of points per cluster
				for(UInt i=0; i < cluster_size.size(); ++i)
				{
					cluster_size[i] = best_n_clusters[i].size();
				}

				//--------------------------------------------------------------
				// fill in cluster_id and cluster_size in SILACData structure
				//--------------------------------------------------------------
				for(UInt i=0; i < best_n_clusters.size(); ++i)
				{
					for(UInt j=0; j < best_n_clusters[i].size(); ++j)
					{
						data[best_n_clusters[i][j]].cluster_id = i;
						data[best_n_clusters[i][j]].cluster_size = best_n_clusters[i].size();
					}
				}
				std::sort(data.begin(),data.end());
				std::reverse(data.begin(),data.end()); // largest clusters first

//~ /*
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
				cluster_size = std::vector<Int>(best_n,0);
				k = -1;
				for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
				{
					++cluster_size[it->cluster_id];
				}
//~ */

				//--------------------------------------------------------------
				//create consensus features
				//--------------------------------------------------------------
				if (out!="")
				{
					for (UInt i=0; i<best_n; ++i)
					{
						DoubleReal rt = 0.0;
						DoubleReal mz = 0.0;
						DoubleReal int_l = 0.0;
						DoubleReal int_h = 0.0;
						// intensity vectors used for linear regression
						std::vector<DoubleReal> i1(3*cluster_size[i]);
						std::vector<DoubleReal> i2(3*cluster_size[i]);
						UInt j=0;
						for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
						{
							if ((UInt)it->cluster_id==i)
							{
								i1[3*j] = it->int1;
								i2[3*j] = it->int4;
								i1[3*j+1] = it->int2;
								i2[3*j+1] = it->int5;
								i1[3*j+2] = it->int3;
								i2[3*j+2] = it->int6;

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
								if (it->int4>int_h)
								{
									int_h = it->int4;
								}
								if (it->int5>int_h)
								{
									int_h = it->int5;
								}
								if (it->int5>int_h)
								{
									int_h = it->int6;
								}
								++j;
							}
						}
						rt /= (DoubleReal)(cluster_size[i]); // average retention time
						Math::LinearRegression linear_reg;
						linear_reg.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i2.begin());
						//create consensus feature
						ConsensusFeature pair;
						pair.setRT(rt);
						pair.setMZ(mz);
						pair.setIntensity(linear_reg.getSlope());
						pair.setCharge(charge);
						pair.setQuality(linear_reg.getRSquared());
						FeatureHandle handle;
						handle.setRT(rt);
						handle.setMZ(mz);
						handle.setIntensity(int_l);
						handle.setCharge(charge);
						handle.setMapIndex(0);
						handle.setElementIndex(i);
						pair.insert(handle);
						handle.setRT(rt);
						handle.setMZ(mz+envelope_distance);
						handle.setIntensity(int_h);
						handle.setCharge(charge);
						handle.setMapIndex(1);
						handle.setElementIndex(i);
						pair.insert(handle);
						all_pairs.push_back(pair);
					}
				}


				//--------------------------------------------------------------
				//create features (for visualization)
				//--------------------------------------------------------------
				if (out_visual!="")
				{
					std::vector<String> colors;
					// greens
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
				}





				//-------------------------------------------------------------
  	    // generate debug output
  	    //-------------------------------------------------------------
				if (getFlag_("silac_debug"))
				{
  	    	String debug_suffix = "_" + String(0.01*floor(mass_separation*100+0.5)) + "Da_" + String(charge) +"+";
  	    	// names of dat files
  	    	String debug_dat = debug_trunk + debug_suffix + ".dat";
  	    	String debug_clusters_dat = debug_trunk + debug_suffix + "_clusters.dat";

  	    	// write all cluster data points to *_clusters.dat
					std::ofstream stream_clusters(debug_clusters_dat.c_str());
					stream_clusters << "cluster_id cluster_size rt mz int1 int2 int3 int4 int5 int6" << std::endl;
					Int current_id = -1;
					for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
					{
						if (it->cluster_id != current_id) stream_clusters << std::endl << std::endl;
						stream_clusters << it->cluster_id << " " << it->cluster_size << " "<< it->rt << " " << it->mz << " " << it->int1 << " " << it->int2 << " " << it->int3 << " " << it->int4 << " " << it->int5 << " " << it->int6 << std::endl;
						current_id = it->cluster_id;
					}
	      	stream_clusters.close();

	      	// write ratios of all cluster to *.dat
					std::ofstream stream_ratios(debug_dat.c_str());
					stream_ratios << "cluster_id cluster_size rt mz ratio" << std::endl;
	      	for (UInt i=0; i<best_n;++i)
					{
						DoubleReal rt = 0.0;
						DoubleReal mz = 0.0;
						DoubleReal int_l = 0.0;
						DoubleReal int_h = 0.0;
						// intensity vectors used for linear regression
						std::vector<DoubleReal> i1(3*cluster_size[i]);
						std::vector<DoubleReal> i2(3*cluster_size[i]);
						UInt j=0;
		   			for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
						{
							if ((UInt)it->cluster_id==i)
							{
								i1[3*j] = it->int1;
								i2[3*j] = it->int4;
								i1[3*j+1] = it->int2;
								i2[3*j+1] = it->int5;
								i1[3*j+2] = it->int3;
								i2[3*j+2] = it->int6;

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
								if (it->int4>int_h)
								{
									int_h = it->int4;
								}
								if (it->int5>int_h)
								{
									int_h = it->int5;
								}
								if (it->int5>int_h)
								{
									int_h = it->int6;
								}
								++j;
							}
						}
						rt /= (DoubleReal)(cluster_size[i]); // average retention time
						Math::LinearRegression linear_reg;
						linear_reg.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i2.begin());
						stream_ratios << i << " " << cluster_size[i] << " " << rt << " " << mz << " " << linear_reg.getSlope() << std::endl;
					}
					stream_ratios.close();
				}

			} //end iterate over all charge states


			//--------------------------------------------------------------
			//Store output
			//--------------------------------------------------------------
  		if (out!="")
  		{
	  		ConsensusXMLFile c_file;
	  		c_file.store(out,all_pairs);
  		}

  		if (out_visual!="")
  		{
  			FeatureXMLFile f_file;
  			f_file.store(out_visual,all_cluster_points);
      }

			//--------------------------------------------------------------
			//write gnuplot script
			//--------------------------------------------------------------
			if (getFlag_("silac_debug"))
			{
				// first lines of the gnuplot script
				String debug_gnuplotscript = debug_trunk + ".input";
				std::ofstream stream_gnuplotscript(debug_gnuplotscript.c_str());
				stream_gnuplotscript << "set terminal postscript eps enhanced colour" << std::endl;
				stream_gnuplotscript << "set size 2.0, 2.0" << std::endl;
				stream_gnuplotscript << "set size square" << std::endl << std::endl;
				//iterate over all charge states
				for (UInt charge=charge_min; charge<=charge_max; ++charge)
				{
  	    	String debug_suffix = "_" + String(0.01*floor(mass_separation*100+0.5)) + "Da_" + String(charge) +"+";
  	    	// names of dat files
  	    	String debug_dat = debug_trunk + debug_suffix + ".dat";
  	    	String debug_clusters_dat = debug_trunk + debug_suffix + "_clusters.dat";
  	    	// names of postscript files
  	    	String debug_ratios = debug_trunk + debug_suffix + "_ratios.eps";
  	    	String debug_sizes = debug_trunk + debug_suffix + "_sizes.eps";
  	    	String debug_clusters = debug_trunk + debug_suffix + "_clusters.eps";
  	    	String debug_Clusters = debug_trunk + debug_suffix + "_Clusters.eps";
  	    	String debug_clustersInt = debug_trunk + debug_suffix + "_clustersInt.eps";
  	    	String debug_ClustersInt = debug_trunk + debug_suffix + "_ClustersInt.eps";

  	    	// write *_clusters.eps
					stream_gnuplotscript << "set output \"" + debug_clusters + "\"" << std::endl;
					stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << ", mass separation = " << String(0.01*floor(mass_separation*100+0.5)) << " Da, charge = " << charge << "+, intensity cutoff = " << intensity_cutoff << ", rt scaling = " + String(0.01*floor(rt_scaling*100+0.5)) + ", cluster number scaling = " + String(0.01*floor(cluster_number_scaling*100+0.5)) + "\"" << std::endl;
					stream_gnuplotscript << "set xlabel \'m/Z (Th)\'" << std::endl;
					stream_gnuplotscript << "set ylabel \'RT (s)\'" << std::endl;
					stream_gnuplotscript << "plot";
		  		for (int i = 0; i < cluster_number[charge]; i++)// iterate over clusters
		  		{
		  			if (i != 0) stream_gnuplotscript << ",";
		  			stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 4:3 title \"cluster " + String(i) + "\"";
		  		}
					stream_gnuplotscript << std::endl;

  	    	// write *_Clusters.eps
					stream_gnuplotscript << "set output \"" + debug_Clusters + "\"" << std::endl;
					stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << ", mass separation = " << String(0.01*floor(mass_separation*100+0.5)) << " Da, charge = " << charge << "+, intensity cutoff = " << intensity_cutoff << ", rt scaling = " + String(0.01*floor(rt_scaling*100+0.5)) + ", cluster number scaling = " + String(0.01*floor(cluster_number_scaling*100+0.5)) + "\"" << std::endl;
					stream_gnuplotscript << "set xlabel \'m/Z (Th)\'" << std::endl;
					stream_gnuplotscript << "set ylabel \'RT (s)\'" << std::endl;
					stream_gnuplotscript << "plot";
		  		for (int i = cluster_min; i <= cluster_max; i++)// iterate over clusters
		  		{
		  			if (i != 0) stream_gnuplotscript << ",";
		  			stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 4:3 title \"cluster " + String(i) + "\"";
		  		}
					stream_gnuplotscript << std::endl;

  	    	// write *_clustersInt.eps
					stream_gnuplotscript << "set output \"" + debug_clustersInt + "\"" << std::endl;
					stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << ", mass separation = " << String(0.01*floor(mass_separation*100+0.5)) << " Da, charge = " << charge << "+, intensity cutoff = " << intensity_cutoff << ", rt scaling = " + String(0.01*floor(rt_scaling*100+0.5)) + ", cluster number scaling = " + String(0.01*floor(cluster_number_scaling*100+0.5)) + "\"" << std::endl;
					stream_gnuplotscript << "set xlabel \'intensity at m/Z\'" << std::endl;
					stream_gnuplotscript << "set ylabel \'intensity at m/Z + " + String(0.01*floor(mass_separation/charge*100+0.5)) + "Th\'" << std::endl;
					stream_gnuplotscript << "plot";
		  		for (int i = 0; i < cluster_number[charge]; i++)// iterate over clusters
		  		{
		  			if (i != 0) stream_gnuplotscript << ",";
		  			stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 5:8 title \"cluster " + String(i) + "\"";
		  		}
					stream_gnuplotscript << std::endl;

  	    	// write *_ClustersInt.eps
					stream_gnuplotscript << "set output \"" + debug_ClustersInt + "\"" << std::endl;
					stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << ", mass separation = " << String(0.01*floor(mass_separation*100+0.5)) << " Da, charge = " << charge << "+, intensity cutoff = " << intensity_cutoff << ", rt scaling = " + String(0.01*floor(rt_scaling*100+0.5)) + ", cluster number scaling = " + String(0.01*floor(cluster_number_scaling*100+0.5)) + "\"" << std::endl;
					stream_gnuplotscript << "set xlabel \'intensity at m/Z\'" << std::endl;
					stream_gnuplotscript << "set ylabel \'intensity at m/Z + " + String(0.01*floor(mass_separation/charge*100+0.5)) + "Th\'" << std::endl;
					stream_gnuplotscript << "plot";
		  		for (int i = cluster_min; i <= cluster_max; i++)// iterate over clusters
		  		{
		  			if (i != 0) stream_gnuplotscript << ",";
		  			stream_gnuplotscript << " \'" + debug_clusters_dat + "\' index " + String(i+1) +" using 5:8 with lines title \"cluster " + String(i) + "\"";
		  		}
					stream_gnuplotscript << std::endl;

  	    	// write *_ratios.eps
					stream_gnuplotscript << "set output \"" + debug_ratios + "\"" << std::endl;
					stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << ", mass separation = " << String(0.01*floor(mass_separation*100+0.5)) << " Da, charge = " << charge << "+, intensity cutoff = " << intensity_cutoff << ", rt scaling = " + String(0.01*floor(rt_scaling*100+0.5)) + ", cluster number scaling = " + String(0.01*floor(cluster_number_scaling*100+0.5)) + "\"" << std::endl;
					stream_gnuplotscript << "set nokey" << std::endl;
					stream_gnuplotscript << "set xlabel \'m/Z\'" << std::endl;
					stream_gnuplotscript << "set ylabel \'ratio\'" << std::endl;
					stream_gnuplotscript << "plot \'" + debug_dat + "\' using 4:5";
					stream_gnuplotscript << std::endl;

  	    	// write *_sizes.eps
					stream_gnuplotscript << "set output \"" + debug_sizes + "\"" << std::endl;
					stream_gnuplotscript << "set title \"SILACAnalyzer " << version_ << ", sample = " << debug_trunk << ", mass separation = " << String(0.01*floor(mass_separation*100+0.5)) << " Da, charge = " << charge << "+, intensity cutoff = " << intensity_cutoff << ", rt scaling = " + String(0.01*floor(rt_scaling*100+0.5)) + ", cluster number scaling = " + String(0.01*floor(cluster_number_scaling*100+0.5)) + "\"" << std::endl;
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
