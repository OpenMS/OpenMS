// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>
#include <OpenMS/SYSTEM/StopWatch.h>

#include <fstream>

using namespace std;

namespace OpenMS
{

  FeaFiTraits::FeaFiTraits() 
  {
  }


  FeaFiTraits::~FeaFiTraits() 
  {	
  }

  void FeaFiTraits::getNextRt(IDX& index) throw (NoSuccessor, Exception::Precondition)
  {
  	//Corrupt index
  	OPENMS_PRECONDITION(index.first<map_.size(), "Scan index outside of map!");
    OPENMS_PRECONDITION(index.second<map_[index.first].size(), "Peak index outside of scan!");
		
		//last scan
    if (index.first == map_.size()-1)
    {
    	throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getNextRt", index);
		}
		
		// perform binary search to find the neighbour in rt dimension
		CoordinateType mz_pos = map_[index.first][index.second].getPos();	// mz value we want to find
		++index.first;
		MapType::SpectrumType::ConstIterator it = lower_bound(map_[index.first].begin(), map_[index.first].end(), map_[index.first-1][index.second], MapType::SpectrumType::PeakType::PositionLess());	
		
		// if the found peak is at the end of the spectrum, there is not much we can do...
		if ( it == map_[index.first].end() )
		{
			// check for empty scans
			if ( map_[index.first].size() > 0 )
	 			index.second = map_[index.first].size()-1;
			else
				index.second = 0;
		}
		// if the found peak is at the beginning of the spectrum, there is also not much we can do ! 
		else if ( it == map_[index.first].begin() ) 
		{
			index.second = 0;
		}
		// see if the next smaller one fits better
		else 
		{	
			// peak to the right is closer (in m/z dimension)
			if (it->getPos() - mz_pos < mz_pos - (it-1)->getPos() )
			{				
				index.second = it - map_[index.first].begin(); 
			}
			else	// left one is closer
			{
				index.second = --it - map_[index.first].begin(); 
			}
		}
  }
  
  
	void FeaFiTraits::getPrevRt(IDX& index) throw (NoSuccessor, Exception::Precondition)
  {
  	//Corrupt index
  	OPENMS_PRECONDITION(index.first<map_.size(), "Scan index outside of map!");
    OPENMS_PRECONDITION(index.second<map_[index.first].size(), "Peak index outside of scan!");
		
		// first scan
		if (index.first == 0)
		{
			throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevRt", index);
		}
		
		// perform binary search to find the neighbour in rt dimension
		CoordinateType mz_pos = map_[index.first][index.second].getPos();
		--index.first;
		MapType::SpectrumType::ConstIterator it = lower_bound(map_[index.first].begin(), map_[index.first].end(), map_[index.first+1][index.second], MapType::SpectrumType::PeakType::PositionLess());	
		
		// if the found peak is at the end of the spectrum, there is not much we can do.
		if ( it == map_[index.first].end() )
		{
	 		// check for empty scans
			if ( map_[index.first].size() > 0 )
	 			index.second = map_[index.first].size()-1;
			else
				index.second = 0;
		}
		// if the found peak is at the beginning of the spectrum, there is not much we can do.
		else if ( it == map_[index.first].begin() ) 
		{
			index.second = 0;
		}
		// see if the next smaller one fits better
		else 
		{	
			// peak to the right is closer (in m/z dimension)
			if (it->getPos() - mz_pos < mz_pos - (it-1)->getPos() )
			{
				index.second = it - map_[index.first].begin(); 
			}
			else
			{
				index.second = --it - map_[index.first].begin(); 
			}
		}
  }
	
	//Calculates the convex hull of a index set and adds it to the feature
	void FeaFiTraits::addConvexHull(const IndexSet& set, DFeature<2>& f) const
	{
		vector< DPosition<2> > points;
		points.reserve(set.size());
		PositionType2D tmp;
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
    {
    	tmp[MZ] = map_[it->first][it->second].getPos();
    	tmp[RT] = map_[it->first].getRetentionTime();
    	points.push_back(tmp);
    }
		f.getConvexHulls().resize(f.getConvexHulls().size()+1);
		f.getConvexHulls()[f.getConvexHulls().size()-1] = points;	
	}

  const DFeatureMap<2>& FeaFiTraits::run(const vector<BaseSeeder*>& seeders, const vector<BaseExtender*>& extenders, const vector<BaseModelFitter*>& fitters)
  {
    // Visualize seeds and extension in TOPPView:
    // get all Seeds and save corresponding peaks as "features"
    // get convex hull of the extension and use it for the "feature"

    // counts the number of features collected so far,
    // is needed for the gnuplot output.
#ifdef DEBUG_FEATUREFINDER
    int nr_feat = 0;
#endif

    if (map_.getSize() == 0)
    {
      cout << "No data provided! Aborting..." << endl;
      return features_;
    }

    if (seeders.size() == 0 || extenders.size() == 0 || fitters.size() == 0)
    {
      cout << "No modules set. Aborting..." << endl;
      return features_;
    }

    // gather information for fitting summary
    map<String,int> exception;									//count exceptions
    int no_exceptions = 0;
    map<String,int> mz_model;									//count used mz models
    map<float,int> mz_stdev;										//count used mz standard deviations
    vector<int> charge(10);											//count used charges
    double corr_mean=0.0, corr_max=0.0, corr_min=1.0; 	//boxplot for correlation

    StopWatch watch;
    unsigned int seed_count = 0;
    try
  	{
      while (true)
      {
				cout << "(1) Seeding ( seed # " << ++seed_count << ")..." << endl;
				IndexSet seed_region = seeders[0]->nextSeed();
				cout << "(2) Extension ..." << endl;
        watch.start();
        IndexSet peaks = extenders[0]->extend(seed_region);
        watch.stop();
        cout << "Time spent for extension: " << watch.getClockTime() << endl;
        watch.reset();
				cout << "(3) ModelFitting ..." << endl;
        try
        {
          watch.start();
          features_.push_back(fitters[0]->fit(peaks));
          watch.stop();
          cout << "Time spent for fitting: " << watch.getClockTime() << endl;
          watch.reset();

#ifdef DEBUG_FEATUREFINDER
          writeGnuPlotFile_(peaks,false,nr_feat++);
#endif
          // gather information for fitting summary
          const DFeature<2>& f = features_.back();

          float corr = f.getOverallQuality();
          corr_mean += corr;
          if (corr<corr_min) corr_min = corr;
          if (corr>corr_max) corr_max = corr;

          // count estimated charge states
          unsigned int ch = f.getCharge();
          if (ch>= charge.size())
          {
          	charge.resize(ch);
          }
        	charge[ch]++;

          const Param& p = f.getModelDescription().getParam();
          ++mz_model[ p.getValue("MZ") ];

          DataValue dp = p.getValue("MZ:isotope:stdev");
          if (dp != DataValue::EMPTY)
          {
          	++mz_stdev[dp];
          }
        }
        catch( BaseModelFitter::UnableToFit ex)
        {
          // set unused flag for all data points
          for (IndexSet::const_iterator it=peaks.begin(); it!=peaks.end(); ++it)
          {
          	getPeakFlag(*it) = FeaFiTraits::UNUSED;
          }
          cout << " " << ex.what() << endl;
          watch.stop();
          cout << "Time spent for fitting: " << watch.getClockTime() << endl;
          watch.reset();
          ++no_exceptions;
          ++exception[ex.getName()];
					
				#ifdef DEBUG_FEATUREFINDER
          writeGnuPlotFile_(peaks,false,nr_feat++);
				#endif
        }

      } // end of while(true)
  	}
    catch(NoSuccessor ex)
    {
    }

    // Print summary:
    Size size = features_.size();

    cout << size << " features were found. " << endl;

//		cout << "seed count " << seed_count << endl;

    cout << "FeatureFinder summary:\n"
    << "Correlation:\n\tminimum: " << corr_min << "\n\tmean: " << corr_mean/size
    << "\n\tmaximum: " << corr_max << endl;

    cout << "Exceptions:\n";
    for (map<String,int>::const_iterator it=exception.begin(); it!=exception.end(); ++it)
    {
      cout << "\t" << it->first << ": " << it->second*100/no_exceptions << "% (" << it->second << ")\n";
    }

    cout << "Chosen mz models:\n";
    for (map<String,int>::const_iterator it=mz_model.begin(); it!=mz_model.end(); ++it)
    {
      cout << "\t" << it->first << ": " << it->second*100/size << "% (" << it->second << ")\n";
    }

    cout << "Chosen mz stdevs:\n";
    for (map<float,int>::const_iterator it=mz_stdev.begin(); it!=mz_stdev.end(); ++it)
    {
      cout << "\t" << it->first << ": " << it->second*100/(size-charge[0]) << "% (" << it->second << ")\n";
    }

    cout << "Charges:\n";
    for (unsigned int i=1; i<charge.size(); ++i)
    {
      if (charge[i]!=0)
      {
        cout << "\t+" << i << ": " << charge[i]*100/(size-charge[0]) << "% (" << charge[i] << ")\n";
      }
		}
#ifdef DEBUG_FEATUREFINDER
    IndexSet empty;
    writeGnuPlotFile_(empty,true,nr_feat);
#endif

    return features_;

  } // end of run(seeders, extenders, fitters)

  void FeaFiTraits::writeGnuPlotFile_(IndexSet peaks, bool last,int nr_feat)
  {
    // write feature + surrounding region to file
    if (!last)
    {
      String gp_fname("plot.gp");
      ofstream gpfile( gp_fname.c_str(), ios_base::app  );

      String file    = String("region") + nr_feat;

			if (nr_feat == 0)
			{
				gpfile << "splot \'" << file << "\' w i title \"\" " << endl;
			}
			else
			{
				gpfile << "replot \'" << file << "\' w i title \"\" " << endl;
			}
      ofstream myfile(file.c_str()); // data file
      IndexSet::const_iterator citer = peaks.begin();

      while (citer != peaks.end())
      {
        myfile << getPeakRt(*citer) << " " << getPeakMz(*citer) << " " << getPeakIntensity(*citer) << endl;
        citer++;
      }
      myfile.close();
      gpfile.close();
    }
    else
    {
      String gp_fname("plot.gp");

      ofstream gpfile( gp_fname.c_str() , ios_base::app );
      gpfile << "pause -1 \'Please hit enter to continue....\' " << endl;
      gpfile.close();
    }
  }


} // end of namespace OpenMS

