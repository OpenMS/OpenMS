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

extern "C"
{
	#include <cluster.h>
}

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
   @page SILACAnalyzer SILACAnalyzer
 
   @brief Determines the ratio of peak pairs in LC-MS data.

   (1) data reduction
   (2) hierarchical clustering in RT-m/Z plane, determine cluster number by maximizing the average silhouette width
   (3) determine intensity ratios by linear regression for each cluster
 
   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

//TODO:
//- several charges
//- documentation
//- debug output and gnuplot scripts


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
        : TOPPBase("SILACAnalyzer","Determination of peak ratios in LC-MS data","0.6.1")
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
			registerDoubleOption_("mz_step_width","<double>",0.01,"step width with which the (interpolated) spectrum is scanned, [m/Z]=Th",false,true);
			setMinFloat_("mz_step_width",0.0);
			registerDoubleOption_("rt_scaling","<double>",0.05,"scaling factor of retention times (Cluster height [s] an\ncluster width [Th] should be of the same order. The clustering algorithms work better for\nsymmetric clusters.)",false,true);
			setMinFloat_("rt_scaling",0.0);
			registerDoubleOption_("cluster_number_scaling","<double>",1.0,"scaling factor of the number of clusters (The average-silhouette-width\nalgorithm returns an 'optimal' number of clusters. This number might need\nto be adjusted by this factor.)",false,true);
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
      DoubleReal cluster_number_scaling = getDoubleOption_("cluster_number_scaling");
			UInt cluster_min = getIntOption_("cluster_min");
			UInt cluster_max = getIntOption_("cluster_max");
      
      String in = getStringOption_("in");
			String out = getStringOption_("out");
			String out_visual = getStringOption_("out_visual");
			
			//-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------
      MzDataFile file;
      MSExperiment<Peak1D> exp;
      
      file.setLogType(log_type_);
      file.load(in,exp);
			
			//output variables
			ConsensusMap output;
			output.getFileDescriptions()[0].filename = in;
			output.getFileDescriptions()[0].label = "light";
			output.getFileDescriptions()[0].size = 0;
			output.getFileDescriptions()[1].filename = in;
			output.getFileDescriptions()[1].label = "heavy";
			output.getFileDescriptions()[1].size = 0;
			FeatureMap<> output_cluster;
      
			
			//iterate over all charges
			for (UInt charge=charge_min; charge<=charge_max; ++charge)
			{
				DoubleReal isotope_distance = 1.0 / (DoubleReal)charge;
				DoubleReal envelope_distance = mass_separation / (DoubleReal)charge;
	      //-------------------------------------------------------------
	      // build SILACData structure
	      //-------------------------------------------------------------		
	      ProgressLogger logger_;
	      std::vector<SILACData> data;
	      
				logger_.setLogType(log_type_);
				logger_.startProgress(0,exp.size(),"reducing raw data");
		    //----------------------------------------------------------------
		    // scan over the entire experiment and write to data structure
		    //----------------------------------------------------------------
		    for (MSExperiment<>::Iterator rt_it=exp.begin(); rt_it!=exp.end(); ++rt_it)
			  {
			  	logger_.setProgress(rt_it-exp.begin());
					Int number_data_points = rt_it->size();
					// read OpenMS data into GSL structure
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
	     		gsl_interp_accel *acc = gsl_interp_accel_alloc();
	     		gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, number_data_points);
	     		gsl_spline_init(spline, &*mz_vec.begin(), &*intensity_vec.begin(), number_data_points);
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
						//bool cond2 = (int_spline3 <= int_spline1) && (int_spline5 <= int_spline3) && (int_spline4 <= int_spline2) && (int_spline6 <= int_spline4); // isotopic peaks within one envelop decrease
						if (cond1)
						{
			  			data.push_back(SILACData(rt_it->getRT(),mz,int_spline1,int_spline3,int_spline5,int_spline2,int_spline4,int_spline6));
						}
					}
				}
				exp.clear();
				logger_.endProgress();
				Int size = data.size(); // number of data points after the reduction
				
	      //-------------------------------------------------------------
	      // generate distance matrix
	      //-------------------------------------------------------------
	      std::vector< std::vector<DoubleReal> > distanceMatrix;
				for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
				{
				  std::vector<DoubleReal> vec;
					for (std::vector<SILACData>::iterator it2=data.begin(); it2!= data.end(); ++it2)
					{
						vec.push_back( sqrt( (it->rt-(*it2).rt)*(it->rt-(*it2).rt)*rt_scaling*rt_scaling+(it->mz-(*it2).mz)*(it->mz-(*it2).mz) ) ); // shrink RT by factor rt_scaling in order to make clusters more symmetric
					}
					distanceMatrix.push_back(vec);
				}
				
	      //-------------------------------------------------------------
	      // copy distance matrix
	      //-------------------------------------------------------------
				double** dm; // distance matrix for the clustering algorithm (will be messed up when tree is generated)
				dm = (double**) malloc(size*sizeof(double*));
			  for (int i = 0; i < size; i++)
			  {
	  			dm[i] = (double*) malloc(size*sizeof(double));
	  		}
	  		for (Int i = 0; i < size; i++)
	  		{
	  			for (Int j = 0; j < size; j++)
	  			{
	  				dm[i][j] = distanceMatrix[i][j];
	  			}
	  		}
				
	      //--------------------------------------------------------------
	      // generate tree
	      //--------------------------------------------------------------
				Node* tree = treecluster(size,size,0,0,0,0,'e','a',dm); // e=EuclideanMetric a=ArithmeticMean=AverageLinkage
	  		int* clusterid = (int*) malloc(size*sizeof(int));
	  		
	      //----------------------------------------------------------------
	      // find number of clusters that maximizes average silhouette width
	      //----------------------------------------------------------------
				logger_.setLogType(log_type_);
				logger_.startProgress(0,size/10,"determining number of clusters");
	  		std::vector<DoubleReal> s(size,0); // average silhoutte width for each number of clusters
	  		for (Int i = 1; i < (Int) size/10; i++) // iterate through (number of clusters)/10 (method not perfect, hence 1/10)
				{
					logger_.setProgress(i);
	  			cuttree(size, tree, i, clusterid);
	  			std::vector<Int> cluster_size(i,0); // number of points per cluster
	  			for (Int j = 0; j < size; j++)// iterate through all points
	  			{
	  				cluster_size[clusterid[j]]++; // count the number of points per cluster
	  			}
	  			for (Int j = 0; j < size; j++)// iterate through all points
		  		{
		  			std::vector<DoubleReal> c(i,0); // average distance of point j to each cluster
	  				for (Int k = 0; k < size; k++)// iterate again through all points
	  				{
	  					c[clusterid[k]] += distanceMatrix[j][k];
	  				}
	  				for (Int k = 0; k < i; k++)
	  				{
	  					c[k] = c[k]/cluster_size[k];
	  				}
	  				DoubleReal a = c[clusterid[j]];
	  				DoubleReal b = 0;
	  				for (Int k = 0; k < i; k++)// iterate over all clusters
	  				{
	  					if (b == 0 || (c[k] < b && k!=clusterid[j])) b = c[k]; // find the nearest one to point j
	  				}
	  				if (cluster_size[clusterid[j]]>1) s[i] += (b-a)/std::max(a,b); // add the silhouette width of point j (if j is the only point in its cluster then silhouette width =0 by definition)
	  			}
	  			s[i] = s[i]/size;
	  		}
				logger_.endProgress();   
	  		
	  		DoubleReal s_max = -1;
	  		Int best_n = 1;
	  		for (Int i = 1; i < (Int) size/10; ++i)
	  		{
	  			if (s[i]>s_max)
	  			{
	  				s_max = s[i];
	  				best_n = i;
	  			}			
	  		}
	  		
	  		//cout << "old: " << best_n << std::endl;
	  		best_n = cluster_number_scaling * best_n; // slightly increase cluster number
	  		//cout << "new: " << best_n << std::endl;
	  		 		
	  		cuttree(size, tree, best_n, clusterid);
	  		
	      //-------------------------------------------------------------
	      // count data points in each cluster
	      //-------------------------------------------------------------  		
				std::vector<Int> cluster_size(best_n,0); // number of points per cluster
				for (Int j = 0; j < size; ++j)// iterate through all points
				{
					cluster_size[clusterid[j]]++; // count the number of points per cluster
				}
	
	      //--------------------------------------------------------------
	      // fill in cluster_id and cluster_size in SILACData structure
	      //--------------------------------------------------------------
	  		Int k = 0;
				for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
				{
					it->cluster_id = clusterid[k];
					it->cluster_size = cluster_size[clusterid[k]];
					++k;
				}
	  		std::sort(data.begin(),data.end());
	  		std::reverse(data.begin(),data.end());
	  		
	      //--------------------------------------------------------------
	      // update cluster_id
	      //--------------------------------------------------------------
				k = -1;
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
				
				//--------------------------------------------------------------
				//create consensus features
				//--------------------------------------------------------------
				if (out!="")
				{
		      for (Int i=0; i<best_n;++i)
					{
						DoubleReal rt = 0.0;
						DoubleReal mz = 0.0;
						DoubleReal int_l = 0.0;
						DoubleReal int_h = 0.0;
						std::vector<DoubleReal> i1(3*cluster_size[i]);
						std::vector<DoubleReal> i2(3*cluster_size[i]);
						UInt j=0;
			   		for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
						{
							if (it->cluster_id==i)
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
						rt /= (DoubleReal)(cluster_size[i]);
						Math::LinearRegression linear_reg;
						linear_reg.computeRegressionNoIntercept(0.95,i1.begin(),i1.end(),i2.begin());
						//create consensus feature
						ConsensusFeature tmp_cluster;
						tmp_cluster.setRT(rt);
						tmp_cluster.setMZ(mz);
						tmp_cluster.setIntensity(linear_reg.getSlope());
						tmp_cluster.setCharge(charge);
						tmp_cluster.setQuality(linear_reg.getRSquared());
						FeatureHandle handle;
						handle.setRT(rt);
						handle.setMZ(mz);
						handle.setIntensity(int_l);
						handle.setCharge(charge);
						handle.setMapIndex(0);
						handle.setElementIndex(i);
						tmp_cluster.insert(handle);
						handle.setRT(rt);
						handle.setMZ(mz+envelope_distance);
						handle.setIntensity(int_h);
						handle.setCharge(charge);
						handle.setMapIndex(1);
						handle.setElementIndex(i);
						tmp_cluster.insert(handle);
						output.push_back(tmp_cluster);
					}
				}
				//--------------------------------------------------------------
				//create features (for visualization)
				//--------------------------------------------------------------
				if (out_visual!="")
				{
					std::vector<String> colors;
					colors.push_back("#000000");
					colors.push_back("#FF0000");
					colors.push_back("#00FF00");
					colors.push_back("#0000FF");
					colors.push_back("#FFFF00");
					colors.push_back("#FF00FF");
					colors.push_back("#00FFFF");

		   		for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
					{
						//light variant
						Feature tmp;
						tmp.setRT(it->rt);
						tmp.setMZ(it->mz);
						DoubleReal intensity = std::max(it->int1,it->int2);
						intensity = std::max(intensity,it->int3);
						tmp.setIntensity(intensity);
						tmp.setCharge(charge);
						tmp.setMetaValue("cluster_id",it->cluster_id);
						tmp.setMetaValue("color",colors[it->cluster_id%colors.size()]);
						output_cluster.push_back(tmp);
					}
				}
				

//	      //--------------------------------------------------------------
//	      // determine file names for debug output
//	      //--------------------------------------------------------------
//				String tmp(in);
//				String out_trunk = tmp.reverse().substr((tmp.suffix('.')).length()+1).reverse() + "_" + String(0.0001*floor(envelope_distance*10000+0.5)) + "Th";
//				String out_clusters_dat = out_trunk + ".clusters.dat";
//				String out_ratios_dat = out_trunk + ".ratios.dat";
//				String out_script = "script.input";
//
//	      //--------------------------------------------------------------
//	      // write SILACData to output stream
//	      //--------------------------------------------------------------
//				std::ofstream stream(out_clusters_dat.c_str());
//				stream << "cluster_id cluster_size rt mz int1 int2 int3 int4 int5 int6" << std::endl;
//				k = -1;
//				for (std::vector<SILACData>::iterator it=data.begin(); it!= data.end(); ++it)
//				{
//					if (it->cluster_id != k) stream << std::endl << std::endl;
//					stream << it->cluster_id << " " << it->cluster_size << " "<< it->rt << " " << it->mz << " " << it->int1 << " " << it->int2 << " " << it->int3 << " " << it->int4 << " " << it->int5 << " " << it->int6 << std::endl;
//					k = it->cluster_id;
//				}
//	      stream.close();
//	      
//	      //-------------------------------------------------------------
//	      // write gnuplot script
//	      //-------------------------------------------------------------  		
//				std::ofstream stream3(out_script.c_str());
//				stream3 << "set terminal postscript eps enhanced colour" << std::endl;
//				stream3 << "set size 2.0,2.0" << std::endl;
//				stream3 << "set size square" << std::endl << std::endl;
//				
//				stream3 << "set output \"" + out_trunk + ".clusters.eps\"" << std::endl;
//				stream3 << "set title \"SILACanalyzer " << version_ << ", intensity_cutoff =" << intensity_cutoff << ", rt_scaling=" + String(0.01*floor(rt_scaling*100+0.5)) + ", cluster_number_scaling=" + String(0.0001*floor(cluster_number_scaling*10000+0.5)) + "\"" << std::endl;
//				stream3 << "set xlabel \'[m/Z]=Th\'" << std::endl;
//				stream3 << "set ylabel \'[RT]=s\'" << std::endl;
//				stream3 << "plot";
//	  		for (int i = 0; i < best_n; i++)// iterate over clusters
//	  		{
//	  			if (i != 0) stream3 << ",";
//	  			stream3 << " \'" + out_clusters_dat + "\' index " + String(i+1) +" using 4:3 title \"cluster " + String(i) + "\"";
//	  		}
//	  		//stream3 << ", \'" + out_trunk + ".manual.dat\' using 3:2 title \"manual\"";
//				
//				stream3 << std::endl << std::endl << "set output \"" + out_trunk + ".clustersInt.eps\"" << std::endl;
//				stream3 << "set title \"SILACanalyzer " << version_ << ", intensity_cutoff =" + String(intensity_cutoff) + ", rt_scaling=" + String(0.01*floor(rt_scaling*100+0.5)) + ", cluster_number_scaling=" + String(0.01*floor(cluster_number_scaling*100+0.5)) + "\"" << std::endl;
//				stream3 << "set xlabel \'intensity at m/Z\'" << std::endl;
//				stream3 << "set ylabel \'intensity at m/Z + " + String(0.0001*floor(envelope_distance*10000+0.5)) + "Th\'" << std::endl;
//				stream3 << "plot";
//	  		for (int i = 0; i < best_n; i++)// iterate over clusters
//	  		{
//	  			if (i != 0) stream3 << ",";
//	  			stream3 << " \'" + out_clusters_dat + "\' index " + String(i+1) +" using 5:8 title \"cluster " + String(i) + "\"";
//	  		}
//				
//				stream3 << std::endl << std::endl << "set output \"" + out_trunk + ".someClusters.eps\"" << std::endl;
//				stream3 << "set title \"SILACanalyzer " << version_ << ", intensity_cutoff =" + String(intensity_cutoff) +  ", rt_scaling=" + String(0.0001*floor(rt_scaling*10000+0.5)) + ", cluster_number_scaling=" + String(0.0001*floor(cluster_number_scaling*10000+0.5)) + "\"" << std::endl;
//				stream3 << "set xlabel \'[m/Z]=Th\'" << std::endl;
//				stream3 << "set ylabel \'[RT]=s\'" << std::endl;
//				stream3 << "plot";
//	  		for (UInt i = cluster_min; i <= cluster_max; i++)// iterate over clusters
//	  		{
//	  			if (i != cluster_min) stream3 << ",";
//	  			stream3 << " \'" + out_clusters_dat + "\' index " + String(i+1) +" using 4:3 title \"cluster " + String(i) + "\"";
//	  		}
//	  		
//				stream3 << std::endl << std::endl << "set output \"" + out_trunk + ".someClustersInt.eps\"" << std::endl;
//				stream3 << "set title \"SILACanalyzer " << version_ << ", intensity_cutoff =" + String(intensity_cutoff) +  ", rt_scaling=" + String(0.0001*floor(rt_scaling*10000+0.5)) + ", cluster_number_scaling=" + String(0.0001*floor(cluster_number_scaling*10000+0.5)) + "\"" << std::endl;
//				stream3 << "set xlabel \'intensity at m/Z\'" << std::endl;
//				stream3 << "set ylabel \'intensity at m/Z + " + String(0.0001*floor(envelope_distance*10000+0.5)) + "Th\'" << std::endl;
//				stream3 << "plot";
//	  		for (UInt i = cluster_min; i <= cluster_max; i++)// iterate over clusters
//	  		{
//	  			if (i != cluster_min) stream3 << ",";
//	  			stream3 << " \'" + out_clusters_dat + "\' index " + String(i+1) +" using 5:8 with lines title \"cluster " + String(i) + "\"";
//	  		}
//	  		
//				stream3 << std::endl << std::endl << "set output \"" + out_trunk + ".ratios.eps\"" << std::endl;
//				stream3 << "set title \"SILACanalyzer " << version_ << ", intensity_cutoff =" + String(intensity_cutoff) +  ", rt_scaling=" + String(0.0001*floor(rt_scaling*10000+0.5)) + ", cluster_number_scaling=" + String(0.0001*floor(cluster_number_scaling*10000+0.5)) + "\"" << std::endl;
//				stream3 << "set xlabel \'m/Z\'" << std::endl;
//				stream3 << "set ylabel \'ratios\'" << std::endl;
//				stream3 << "plot \'" + out_ratios_dat + "\' using 4:5";
//	  		
//				stream3 << std::endl << std::endl << "set output \"" + out_trunk + ".size.eps\"" << std::endl;
//				stream3 << "set title \"SILACanalyzer " << version_ << ", intensity_cutoff =" + String(intensity_cutoff) +  ", rt_scaling=" + String(0.0001*floor(rt_scaling*10000+0.5)) + ", cluster_number_scaling=" + String(0.0001*floor(cluster_number_scaling*10000+0.5)) + "\"" << std::endl;
//				stream3 << "set xlabel \'cluster ID\'" << std::endl;
//				stream3 << "set ylabel \'size\'" << std::endl;
//				stream3 << "plot \'" + out_ratios_dat + "\' using 1:2";
//	  		
//	      stream3.close();
//				
//	  		std::cout << "number of reduced data points: " << size << std::endl;
//	  		std::cout << "number of clusters: " << best_n << std::endl << std::endl;
  		}
  		
			//--------------------------------------------------------------
			//Store output
			//--------------------------------------------------------------
  		if (out!="")
  		{
	  		ConsensusXMLFile c_file;
	  		c_file.store(out,output);
  		}
  		
  		if (out_visual!="")
  		{
  			FeatureXMLFile f_file;
  			f_file.store(out_visual,output_cluster);
      }

			if (getFlag_("silac_debug"))
			{
				std::cout << "Writing debug info" << std::endl;
			}
      
      return EXECUTION_OK;
    }
};

int main(int argc, const char** argv )
{
  TOPPSILACAnalyzer tool;
  return tool.main(argc,argv);
}
