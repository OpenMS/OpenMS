#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FeatureFinderAlgorithmWatershed_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FeatureFinderAlgorithmWatershed_H

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
    @nWatersheds in digital spaces: an efficient algorithm based onimmersion simulations
    @nVincent, L.   Soille, P.  
		@nIEEE Transactions on Pattern Analysis and Machine Intelligence, 1991, 13 (6), 583-598
    
    @experimental Currently only the watershed segmentation is returnes, not real features!
    
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
        Real intensity; //TODO: Is a real intensity good for the runtime? Would a integer range e.g. [0:10000] be better?
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
      
      //debug image
      QImage image_;
      UInt peaks_;
      DoubleReal mz_sampling_;

      std::vector< std::vector<GridPoint> > data_;
      std::vector<GridPoint*> data_ptrs_;
      std::deque<GridPoint*> fifo_;
      	
      bool debug_;

    public: 
      FeatureFinderAlgorithmWatershed() :
        FeatureFinderAlgorithm<PeakType,FeatureType>()
      {
      	//Algorithm parameters
        this->setName("Watershed");
        this->defaults_.setValue("mz_sampling",1.0,"Sampling rate for m/z dimension.");
				this->defaults_.setMinFloat("mz_sampling",0.0);
				//Debug flags
        this->defaults_.setValue("debug","true","run in debug mode");
        this->defaults_.setValidStrings("debug",StringList::create("true,false"));
				this->defaults_.setValue("debug:image_name","debug.png","Image of the resampled input data and the watershed segmentation");
				this->defaults_.setValue("debug:gradient","0,#FFFFFF;35,#888888;100,#000000","Gradient of greyscale image");
        this->defaultsToParam_();
      }
      virtual void run()
      {
        //Label constants
        const Int FICTIOUS = -3;
        const Int MASK = -2;
        const Int INIT = -1;
        const Int WATERSHED = 0;
        
        //---------------------------------------------------------------------------
        //Step 1:
        //Initialisation (every peak gets the init value)
        //---------------------------------------------------------------------------
        debug_ = param_.getValue("debug").toBool();
        mz_sampling_ = (DoubleReal)(param_.getValue("mz_sampling"));
        peaks_ = (map_->getMaxMZ() - map_->getMinMZ()) / mz_sampling_;
        
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
          
          for (UInt p=0; p<(*map_)[s].size(); ++p)
          {
            lip.addValue((*map_)[s][p].getMZ(), (*map_)[s][p].getIntensity());
          }
          std::vector<GridPoint> spectrum_points;
          spectrum_points.reserve(peaks_);
          for (UInt p=0; p<peaks_; ++p )
		    	{
           	GridPoint current_point;
            current_point.spectrum = s;
            current_point.peak = p;
            current_point.intensity = lip.getData().at(p);
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
					
        	//convert the map to a grayscale image (png-format)
			    MultiGradient gradient;
				  gradient.fromString(String("Linear|") + param_.getValue("debug:gradient").toString());
          image_ = QImage(peaks_, data_.size(), QImage::Format_RGB32);;
        	for (UInt s=0;s<data_.size();++s)
        	{
        		for (UInt p=0;p<data_[s].size();++p)
	        	{
              image_.setPixel(p, data_.size() - s - 1, gradient.interpolatedColorAt(100.0 * data_[s][p].intensity/max_int).rgb());
            }
          }
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
          fict.flag = FICTIOUS;
          fifo_.push_back(&fict);        
          UInt current_dist = 1;
          while (true)
          {
            GridPoint* current_point = fifo_.front();
            fifo_.pop_front();
            if (current_point->flag == FICTIOUS)
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

        if (debug_)
        {
					//debug info
        	std::cout << "Labels: " << current_label << std::endl;
        	//paint watersheds on the input image
          for (UInt i=0; i<data_.size(); ++i)
          {
            for (UInt j=0; j<peaks_; ++j)
            {
              if (data_[i][j].flag == WATERSHED)
              {
              	image_.setPixel(j,data_.size() - 1 - i,QColor("RED").rgb());
              }
            }
          }
          //store image
          image_.save(param_.getValue("debug:image_name").toString().toQString(), "PNG");      
        }  
        ff_->endProgress();


        //---------------------------------------------------------------------------
        //Step 3:
        //Create features
        //---------------------------------------------------------------------------
				ff_->startProgress(0, data_ptrs_.size(), "Creating features");
				features_->resize(current_label);
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
        		Feature& feature = (*features_)[point.flag-1];
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
				//calculate convex hulls and remove features without points
				for (Int i=features_->size()-1; i>=0 ; --i)
				{
					if (points[i].empty())
					{
						features_->erase(features_->begin()+i);
					}
					else
					{
						(*features_)[i].getConvexHulls().push_back(points[i]);
						(*features_)[i].setMetaValue("label",i);
						(*features_)[i].setMetaValue("contained_points",(UInt)(points[i].size()));
					}
				}

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
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FeatureFinderAlgorithmWatershed_H
