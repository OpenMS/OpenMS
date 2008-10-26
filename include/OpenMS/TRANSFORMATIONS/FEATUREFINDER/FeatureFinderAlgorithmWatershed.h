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
// $Maintainer: Dominik Damerow $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMWATERSHED_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMWATERSHED_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>
#include <iostream>
#include <fstream>
#include <deque>
#include <math.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>

#include <QtGui/QImage>
#include <QtGui/QColor>

namespace OpenMS
{
  /** 
		@brief FeatureFinderAlgorithm implementation using the Watershed Segmentation.
    
    The watershed segmentation algorithm is based on the paper:
    @n Watersheds in digital spaces: an efficient algorithm based onimmersion simulations
    @n Vincent, L.   Soille, P.  
		@n IEEE Transactions on Pattern Analysis and Machine Intelligence, 1991, 13 (6), 583-598
    
    @experimental Currently only the watershed segmentation is returned, not real features!
    
    @ingroup FeatureFinder
	*/
  template<class PeakType, class FeatureType> class FeatureFinderAlgorithmWatershed
  	: public FeatureFinderAlgorithm<PeakType, FeatureType>,
    	public FeatureFinderDefs
  {
    protected:
 
      using FeatureFinderAlgorithm<PeakType, FeatureType>::param_;
      using FeatureFinderAlgorithm<PeakType, FeatureType>::map_;
      using FeatureFinderAlgorithm<PeakType, FeatureType>::features_;
      using FeatureFinderAlgorithm<PeakType, FeatureType>::ff_;
    
			/// Internal respresentation of (resampled) data points
      struct GridPoint
      {
        UInt spectrum;
        UInt peak;
        UInt intensity;
        UInt distance;
        Int flag;
        
      	///Comparison functor needed to sort GridPoints
	      struct GridPointLess
	      	: std::binary_function<GridPoint,GridPoint,bool>
	      {
	      	bool operator()(const GridPoint& a,const GridPoint& b) const
					{
						return a.intensity<b.intensity;
					}
	      };
      };
      
      UInt peaks_;
      DoubleReal mz_sampling_;
      DoubleReal cutoff_factor_;

      std::vector< std::vector<GridPoint> > data_;
      std::vector<GridPoint*> data_ptrs_;
      std::deque<GridPoint*> fifo_;
      	
      bool debug_;
      bool apply_log_;

    public: 
      FeatureFinderAlgorithmWatershed() :
        FeatureFinderAlgorithm<PeakType,FeatureType>()
      {
      	//Algorithm parameters
        this->setName("Watershed");
        this->defaults_.setValue("mz_sampling",1.0,"Sampling rate for m/z dimension.");
				this->defaults_.setMinFloat("mz_sampling",0.0);
        this->defaults_.setValue("cutoff_factor",1.0,"Only features with a size of average size/cutoff_factor are allowed.");
        this->defaults_.setMinFloat("cutoff_factor",0.0);
        this->defaults_.setMaxFloat("cutoff_factor",7.0);
        this->defaults_.setValue("apply_log","false","Apply log transformation");
        this->defaults_.setValidStrings("apply_log",StringList::create("true,false"));
				//Debug flags
        this->defaults_.setValue("debug","false","run in debug mode");
        this->defaults_.setValidStrings("debug",StringList::create("true,false"));
        this->defaultsToParam_();
      }
      virtual void run()
      {
        features_->clear();
        
        //Label constants
        const Int FICTITIOUS = -3;
        const Int MASK = -2;
        const Int INIT = -1;
        const Int WATERSHED = 0;
        
        //---------------------------------------------------------------------------
        //Step 1:
        //Initialisation (every peak gets the init value)
        //---------------------------------------------------------------------------
        debug_ = param_.getValue("debug").toBool();
        mz_sampling_ = (DoubleReal)(param_.getValue("mz_sampling"));
        cutoff_factor_ = (DoubleReal)(param_.getValue("cutoff_factor"));
        apply_log_ = param_.getValue("apply_log").toBool();
        peaks_ = UInt( ceil((map_->getMaxMZ() - map_->getMinMZ()) / mz_sampling_));
        Real normalizing_factor = 0;
        if (apply_log_)
        {
          normalizing_factor = 10000/(log(map_->getMaxInt())+1);
        }
        else
        {
          normalizing_factor = 10000/(map_->getMaxInt());
        }
        //---------------------------------------------------------------------------
      	//RESAMPLE AND BUILD MAIN DATASTRUCTURE
				//do linear resampling in m/z dimension
        ff_->startProgress(0, map_->size(), "Resampling of input data");
        data_.reserve(map_->size());
        for (UInt s=0; s<map_->size(); ++s)
        {
        	ff_->setProgress(s);
          if ((*map_)[s].getMSLevel()!=1)
      		{
    				continue;
    			}
    			Math::LinearInterpolation<DoubleReal,DoubleReal> lip;
		      lip.getData().resize(peaks_);
          lip.setMapping( 0, map_->getMinMZ(), peaks_-1, map_->getMaxMZ() );
          
          if (apply_log_)
          {
            for (UInt p=0; p<(*map_)[s].size(); ++p)
            {
              lip.addValue((*map_)[s][p].getMZ(),log(1 + (*map_)[s][p].getIntensity()));
            }
          }
          else
          {
            for (UInt p=0; p<(*map_)[s].size(); ++p)
            {
              lip.addValue((*map_)[s][p].getMZ(),((*map_)[s][p].getIntensity()));
            }
          }
          
          std::vector<GridPoint> spectrum_points;
          spectrum_points.reserve(peaks_);
          for (UInt p=0; p<peaks_; ++p )
		    	{
           	GridPoint current_point;
            current_point.spectrum = s;
            current_point.peak = p;
            current_point.intensity = UInt(round(lip.getData().at(p)*normalizing_factor));
            current_point.distance = 0;
            current_point.flag = INIT;
            spectrum_points.push_back(current_point);
	    		}
	    		data_.push_back(spectrum_points);
        }
        ff_->endProgress();
        
    		//---------------------------------------------------------------------------
        //debug output => png image
        if (debug_)
        {
        	std::cout << "Spectra/peaks: " << data_.size() << "/" << peaks_ << " (overall points: " << data_.size()*peaks_ << ")" << std::endl;
					//determine maximum intensity of resampled data
        	DoubleReal max_int = 0.0;
        	DoubleReal min_int = data_[0][0].intensity;
        	for (UInt s=0;s<data_.size();++s)
        	{
        		for (UInt p=0;p<data_[s].size();++p)
	        	{
	        		if (data_[s][p].intensity > max_int)
	        		{
	        			max_int = data_[s][p].intensity;
	        		}
	        		if (data_[s][p].intensity < min_int)
	        		{
	        			min_int = data_[s][p].intensity;
	        		}
	        	}
        	}
        	std::cout << "min/max intensity: " << min_int << "/" << max_int << std::endl;
        }
        //---------------------------------------------------------------------------
				//SORTING
        //sort seed vector by intensity of peaks (highest first)
       	data_ptrs_.reserve(data_.size()*peaks_);
      	for (UInt s=0;s<data_.size();++s)
      	{
      		for (UInt p=0;p<data_[s].size();++p)
        	{
        		data_ptrs_.push_back(&(data_[s][p]));
          }
        }
        sort(data_ptrs_.begin(),data_ptrs_.end(),reverseComparator(pointerComparator(typename GridPoint::GridPointLess())));

       
        //---------------------------------------------------------------------------
        //Step 2:
        //Flooding step
        //---------------------------------------------------------------------------

        ff_->startProgress(0, data_ptrs_.size(), "Watershed segmentation");

        //labels for basins (>0)
        UInt current_label = 0;
     		
				//Loop over intensities
				UInt i=0;
        while (i<data_ptrs_.size())
        {
          ff_->setProgress(i);
	        //----------------------------------------------------------------------------------
					//LOOK AT ALL POINTS OF THE CURRENT INTENSITY AND PUT THEM IN THE QUEUE IF NECESSARY
          Real current_intensity = data_ptrs_[i]->intensity;          
          //std::cout << "i: " << i << "/" << data_ptrs_.size() << " (intensity: " << current_intensity << ")" << std::endl;
          UInt j = i;
          while (j<data_ptrs_.size() && data_ptrs_[j]->intensity >= current_intensity)
          {
            GridPoint* current_point = data_ptrs_[j];
            if (current_point->intensity == current_intensity)
            {
              current_point->flag = MASK;
            }
            
            std::vector<GridPoint*> neighbors; 
            getNeighbors_(current_point,neighbors);
            for (UInt n=0; n<neighbors.size(); ++n)
            {
              if ( neighbors[n]->flag == WATERSHED || neighbors[n]->flag > 0)
              {
                neighbors[n]->distance = 1;
                fifo_.push_back(neighbors[n]);
              }
            }
            ++j;
          }// end while
          //std::cout << "Points of similar intensity: " << j-i << std::endl;

	        //----------------------------------------------------------------------------------
					//PROCESS THE POINTS IN THE QUEUE
          GridPoint fict;
          fict.flag = FICTITIOUS;
          fifo_.push_back(&fict);        
          UInt current_dist = 1;
          while (true)
          {
            GridPoint* current_point = fifo_.front();
            fifo_.pop_front();
            if (current_point->flag == FICTITIOUS)
            {
              if (fifo_.empty())  break;
              fifo_.push_back(&fict);
              current_dist++;
              current_point = fifo_.front();
              fifo_.pop_front();
            }
            
            // for each labeled or watershed neighbor with a distance < current_dist
            // put neighbors in vector
            std::vector<GridPoint*> neighbors;
            getNeighbors_(current_point,neighbors);
            for (UInt n = 0; n < neighbors.size(); ++n)
            {
              GridPoint* neighbor = neighbors[n];
              if (((neighbor->flag == WATERSHED) || (neighbor->flag > 0)) && (neighbor->distance < current_dist))
              {
                if (neighbor->flag > 0)
                {
                  if (current_point->flag == MASK || current_point->flag == WATERSHED)
                  {
                    current_point->flag = neighbor->flag;
                  }
                  else if (current_point->flag != neighbor->flag)
                  {
                    current_point->flag = WATERSHED;
                  }
                }
                else if (current_point->flag == MASK)
                {
                  current_point->flag = WATERSHED; 
                }          
              }
              else if (neighbor->flag == MASK && neighbor->distance == 0)
              {
                neighbor->distance = current_dist + 1;
                fifo_.push_back(neighbor);
              }
            }//for
          }//while(true)
          
          //----------------------------------------------------------------------------------
          //CHECK IF NEW MINIMA HAVE BEEN DISCOVERED
          for (UInt j2=i; j2<j; ++j2)
          {
            GridPoint* current_point = data_ptrs_[j2];
            //distance is reset to 0 
            current_point->distance = 0;
            if (current_point->flag == MASK)
            {
              fifo_.push_back(current_point);
              //label start at 1 (current_label was initiated with 0)
              ++current_label;
              current_point->flag = current_label;
              while (!fifo_.empty())
              {
                std::vector<GridPoint*> neighbors; 
                getNeighbors_(fifo_.front(),neighbors);
                fifo_.pop_front();
                for (UInt n=0; n<neighbors.size(); ++n)
                {
                  if(neighbors[n]->flag == MASK)
                  {
                    fifo_.push_back(neighbors[n]);
                    neighbors[n]->flag = current_label;
                  }
                }
              }
            }
          }
          
          //move i to j (go to next intensity level)
          i=j;
 				}
        ff_->endProgress();
        
        if (debug_)
        {
					//debug info
        	std::cout << "Labels: " << current_label << std::endl;      
        }  

        //---------------------------------------------------------------------------
        //Step 3:
        //Create features
        //---------------------------------------------------------------------------
				ff_->startProgress(0, data_ptrs_.size(), "Creating features");
				FeatureMap<> tmp_features;
				tmp_features.resize(current_label);
        std::vector<ConvexHull2D::PointArrayType> points;
        points.resize(current_label);
        for (UInt i=0; i<data_ptrs_.size(); ++i)
        {
        	ff_->setProgress(i);
        	const GridPoint& point = *(data_ptrs_[i]);
        	if (point.flag>0)
        	{
        		//calculate RT and m/z position
        		DoubleReal rt = (*map_)[point.spectrum].getRT();
        		DoubleReal mz = map_->getMinMZ() + (0.5+point.peak)*mz_sampling_;
        		//update feature center (to maximum)
        		Feature& feature = tmp_features[point.flag-1];
        		if (point.intensity > feature.getIntensity())
        		{
        			feature.setIntensity(point.intensity);
        			feature.setRT(rt);
        			feature.setMZ(mz);
        		}
        		//add point to convex hull points
        		points[point.flag-1].push_back(ConvexHull2D::PointType(rt,mz));
        	}
        }
        ff_->endProgress();
        
				//Calculate the average contained points
        ff_->startProgress(0, tmp_features.size(),"Calculating average contained points");
        Real average_points = 0;
        UInt counter = 0; 
        UInt pointssize[500];
        for (UInt l = 0; l<500;l++)          
          {  pointssize[l] = 0;
          }       
        for (UInt i=0; i<tmp_features.size(); ++i)
				{
          ff_->setProgress(i);
          if (!points[i].empty())
          {
            counter++;
            average_points += points[i].size();
            if(points[i].size() < 5000)
            {
              pointssize[(UInt)(points[i].size()/10)] +=1;
            }
          }
        }
        average_points /= counter * cutoff_factor_;
         
        ff_->endProgress();
        //calculate convex hulls and copy features with points to the output array
				ff_->startProgress(0, tmp_features.size(), "Calculating feature convex hulls");
				features_->reserve(tmp_features.size()); 
				for (UInt i=0; i<tmp_features.size(); ++i)
				{
					ff_->setProgress(i);
          
					if (points[i].size() > average_points)
					{
						features_->push_back(tmp_features[i]);
						features_->back().getConvexHulls().push_back(points[i]);
						features_->back().setMetaValue("label",i);
						features_->back().setMetaValue("contained_points",(UInt)points[i].size());
            
					}
				}
        std::cout << std::endl;
        ff_->endProgress();
        
        
				//debug info
        if (debug_)
        {
        	std::cout << "Features: " << features_->size() << std::endl;
        }
        
				//---------------------------------------------------------------------------
        //Step 4:
        //Cleaning up
        //---------------------------------------------------------------------------
				data_.clear();
      	data_ptrs_.clear();
      	fifo_.clear();
      }
      
      static FeatureFinderAlgorithm<PeakType,FeatureType>* create()
      {
        return new FeatureFinderAlgorithmWatershed();
      }

      static const String getProductName()
      {
        return "watershed";
      }
		
		protected:
			
			inline void getNeighbors_(GridPoint* point, std::vector<GridPoint*>& neighbors)
			{
				if (point->spectrum != 0)
        {
          neighbors.push_back(&data_[point->spectrum - 1][point->peak]);
        }
        if (point->spectrum < data_.size() - 1)
        {
          neighbors.push_back(&data_[point->spectrum + 1][point->peak]);
        }
        if (point->peak != 0)
        {
          neighbors.push_back(&data_[point->spectrum][point->peak-1]);
        }
        if (point->peak < data_[point->spectrum].size()-1)
        {
          neighbors.push_back(&data_[point->spectrum][point->peak+1]);
        } 
			}
		
    private:
      
      // not implemented
      FeatureFinderAlgorithmWatershed& operator=(const FeatureFinderAlgorithmWatershed&);
      // not implemented      
      FeatureFinderAlgorithmWatershed(const FeatureFinderAlgorithmWatershed&);
  };// FeatureFinderAlgorithmWatershed
}// namespace OpenMS
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMWATERSHED_H
