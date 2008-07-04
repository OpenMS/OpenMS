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
// $Maintainer: Marcel Grunert $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMWAVELET_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMWAVELET_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelFitter.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
	/** 
    @brief FeatureFinderAlgorithm implementation using the IsotopeWavelet and the ModelFitter.

    IsotopeWavelet (Seeding & Extension) and ModelFitter (using EMG in RT dimension and improved IsotopeModel in dimension of mz)

    @ref FeatureFinderAlgorithmWavelet_Parameters are explained on a separate page.
	
    @ingroup FeatureFinder
   */
    template<class PeakType, class FeatureType> 
    class FeatureFinderAlgorithmWavelet :
      public FeatureFinderAlgorithm<PeakType, FeatureType>,
      public FeatureFinderDefs
      {
          public:	  	
            
            ///@name Type definitions
            //@{
            typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::MapType MapType;
            typedef typename MapType::SpectrumType SpectrumType;
            typedef typename PeakType::CoordinateType CoordinateType;
            typedef typename PeakType::IntensityType IntensityType;
            //@}
            
            /// default constructor 
            FeatureFinderAlgorithmWavelet() : FeatureFinderAlgorithm<PeakType,FeatureType>()
            {
							this->defaults_ = getDefaultParameters();
              this->check_defaults_ = false;
              this->defaultsToParam_();
            }
            
            virtual Param getDefaultParameters() const
            {
              Param tmp;
    
              tmp.setValue("max_charge", 2, "The maximal charge state to be considered.", false);
              tmp.setValue("intensity_threshold", 1, "The final threshold t' is build upon the formula: t' = av+t*sd where t is the intensity_threshold, av the average intensity within the wavelet transformed signal and sd the standard deviation of the transform. If you set intensity_threshold=-1, t' will be zero. For single scan analysis (e.g. MALDI peptide fingerprints) you should start with an intensity_threshold around 0..1 and increase it if necessary.", false);
              tmp.setValue("rt_votes_cutoff", 5, "A parameter of the sweep line algorithm. It" "subsequent scans a pattern must occur to be considered as a feature.", false);
              tmp.setValue("rt_interleave", 2, "A parameter of the sweep line algorithm. It determines the maximum number of scans (w.r.t. rt_votes_cutoff) where a pattern is missing.", false);
              tmp.setValue("recording_mode", 1, "Determines if the spectra have been recorded in positive ion (1) or negative ion (-1) mode.", true);
              tmp.setValue("charge_threshold", 0.1, "All features/seeds (found by isotope wavelet) get a set of possible charges. Every charge holds a score and the charge threshold limits the number of charge states to be considered (in ModelFitter).", true);

              ModelFitter<PeakType,FeatureType> fitter(this->map_, this->features_, this->ff_);
              tmp.insert("fitter:", fitter.getParameters());
              tmp.setSectionDescription("fitter", "Settings for the modefitter (Fits a model to the data determinging the probapility that they represent a feature.)");
              
              return tmp;
            }
  
          virtual void run()
          {
            
            ModelFitter<PeakType,FeatureType> fitter(this->map_, this->features_, this->ff_);
            Param params;
            params.setDefaults(this->getParameters().copy("fitter:",true));
            params.setValue("fit_algorithm", "wavelet");
            fitter.setParameters(params);
        
            /// Summary of fitting results
            struct Summary
            {
              std::map<String,UInt> exception; //count exceptions
              UInt no_exceptions;
              std::map<String,UInt> mz_model; //count used mz models
              std::map<float,UInt> mz_stdev; //count used mz standard deviations
              std::vector<UInt> charge; //count used charges
              CoordinateType corr_mean, corr_max, corr_min; 	//boxplot for correlation
              
              /// Initial values
              Summary() :
                  no_exceptions(0),
              corr_mean(0),
              corr_max(0),
              corr_min(1)
              {}
            
            } summary;
  
            //---------------------------------------------------------------------------
            //Step 1:
            //Find seeds with IsotopeWavelet
            //---------------------------------------------------------------------------
            
            CoordinateType max_mz = this->map_->getMax()[1];
            CoordinateType min_mz = this->map_->getMin()[1];
            
            IsotopeWaveletTransform<PeakType> iwt (min_mz, max_mz, max_charge_);
            
            this->ff_->setLogType (ProgressLogger::CMD);
            this->ff_->startProgress (0, 3*this->map_->size(), "analyzing spectra");  

            UInt RT_votes_cutoff = RT_votes_cutoff_;
            //Check for useless RT_votes_cutoff_ parameter
            if (RT_votes_cutoff_ > this->map_->size())
              RT_votes_cutoff = 0;
        
            for (UInt i=0, j=0; i<this->map_->size(); ++i)
            {	
              std::vector<MSSpectrum<PeakType> > pwts (max_charge_, this->map_->at(i));
              std::cout << "Spectrum " << i << " (" << this->map_->at(i).getRT() << ") of " << this->map_->size()-1 << " ... " ; 
              std::cout.flush();
          
              iwt.getTransforms (this->map_->at(i), pwts, max_charge_, mode_);
              this->ff_->setProgress (++j);
    
              #ifdef OPENMS_DEBUG
                std::cout << "transform ok ... "; std::cout.flush();
              #endif
          
              iwt.identifyCharges (pwts,  this->map_->at(i), i, ampl_cutoff_);
              this->ff_->setProgress (++j);
                        
              #ifdef OPENMS_DEBUG
              std::cout << "charge recognition ok ... "; std::cout.flush();
              #endif
              
              iwt.updateBoxStates(*this->map_, i, RT_interleave_, RT_votes_cutoff);
              this->ff_->setProgress (++j);
                      
              #ifdef OPENMS_DEBUG
              std::cout << "updated box states." << std::endl;
              #endif

              std::cout.flush();
            };

            this->ff_->endProgress();
        
            //Forces to empty OpenBoxes_ and to synchronize ClosedBoxes_ 
            iwt.updateBoxStates(*this->map_, INT_MAX, RT_interleave_, RT_votes_cutoff); 
  
            #ifdef OPENMS_DEBUG
            std::cout << "Final mapping."; std::cout.flush();
            #endif
              
          //---------------------------------------------------------------------------
          //Step 2:
          //Calculate bounding box
          //---------------------------------------------------------------------------
          
          // get the closed boxes from IsotopeWavelet
          std::multimap<CoordinateType, Box> boxes = iwt.getClosedBoxes();
 
          // total number of features
          UInt counter_feature = 1;
            
          typename std::multimap<CoordinateType, Box>::iterator iter;
          typename Box::iterator box_iter;
          UInt best_charge_index; CoordinateType c_mz;
          UInt c_charge; // UInt peak_cutoff;
          CoordinateType av_intens=0, av_mz=0;// begin_mz=0; 
          
        	this->ff_->setLogType (ProgressLogger::CMD);
          this->ff_->startProgress (0, boxes.size(), "model fitting ...");  

        	UInt seeds = 0;
        	
          // for all seeds ... 
          for (iter=boxes.begin(); iter!=boxes.end(); ++iter)
          {		
            this->ff_->setProgress (++seeds);
   
            Box& c_box = iter->second;
            std::vector<CoordinateType> charge_votes (max_charge_, 0), charge_binary_votes (max_charge_, 0);
  
  		      //Let's first determine the charge
            for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
            {
              charge_votes[box_iter->second.c] += box_iter->second.score;
              ++charge_binary_votes[box_iter->second.c];
            }
              					
						// Charge voting
            CoordinateType votes = 0;
            for (UInt i=0; i<max_charge_; ++i) votes += charge_votes[i];
            	
            UInt first_charge = 1;
            UInt last_charge = 1; 
            bool set_first = false;
            	
            // get score in percent and set charges
            for (UInt i=0; i<max_charge_; ++i)
            {
              CoordinateType perc_score = charge_votes[i]/votes;
              if (perc_score >= charge_threshold_) 
              {
                if (!set_first) 
                {
                  first_charge = i+1;
                  last_charge = i+1;
                  set_first = true;
                }
            			
                if (last_charge < (i+1)) last_charge = i+1;
              }
            }
            
            // Feature with best correlation
            Feature final_feature;
            
            // quality, correlation for several charges
            CoordinateType quality_feature = 0.0;
            CoordinateType max_quality_feature = -1.0;
      
						//---------------------------------------------------------------------------
         		// Now, check different charges ... 
            //---------------------------------------------------------------------------
            if (first_charge<=last_charge && first_charge >0 && last_charge>0)
            for (UInt i=first_charge; i<=last_charge; ++i)
            {            
              best_charge_index = i-1;
              
            	// Pattern found in too few RT scan 
            	if (charge_binary_votes[best_charge_index] < RT_votes_cutoff && RT_votes_cutoff <= this->map_->size()) 
           	 	{
            		continue;
            	}
                
              // that's the finally predicted charge state for the pattern
              c_charge = best_charge_index + 1; 
              
              //---------------------------------------------------------------------------
              // Get the boundaries for the box with specific charge
              //---------------------------------------------------------------------------
              av_intens=0, av_mz=0;

              // Index set for seed region
              ChargedIndexSet region;
              for (box_iter=c_box.begin(); box_iter!=c_box.end(); ++box_iter)
              {
                c_mz = box_iter->second.mz;
                
                // begin/end of peaks in spectrum
								// peak_cutoff = iwt.getPeakCutOff (c_mz, c_charge);
                //begin_mz = c_mz - NEUTRON_MASS/(CoordinateType)c_charge;
                //const SpectrumType& spectrum = this->map_->at(box_iter->second.RT_index);
								
                UInt spec_index_begin = box_iter->second.MZ_begin; //spectrum.findNearest(begin_mz);
				        UInt spec_index_end = box_iter->second.MZ_end; //spectrum.findNearest(end_mz);
                
                if (spec_index_end >= this->map_->at(box_iter->second.RT_index).size() )
                	break;
                
                // compute index set for seed region
                for (UInt p=spec_index_begin; p<=spec_index_end; ++p)
                {
                  region.insert(std::make_pair(box_iter->second.RT_index,p));
                } 
                
								if (best_charge_index == box_iter->second.c)
                {				
                  av_intens += box_iter->second.intens;
                  av_mz += c_mz*box_iter->second.intens;
                };
                
              };
       
              // calculate the average intensity
              av_intens /= (CoordinateType)charge_binary_votes[best_charge_index];
              // calculate monoisotopic peak
              av_mz /= av_intens*(CoordinateType)charge_binary_votes[best_charge_index];
              // Set charge for seed region
              region.charge_ = c_charge;
              
              //---------------------------------------------------------------------------
              // Step 3:
              // Model fitting
              //---------------------------------------------------------------------------
              try
              {
                // set monoisotopic mz
                fitter.setMonoIsotopicMass(av_mz);
        
								// model fitting
                Feature feature = fitter.fit(region);
                
                // quality, correlation
                quality_feature = feature.getOverallQuality();
  
                if (quality_feature > max_quality_feature)
                {
                  max_quality_feature = quality_feature;
                  final_feature = feature;
                }
                
                // Now, lets see what is the best charge and hence feature
                if (i==last_charge)
                {
                    this->features_->push_back(final_feature);

                    // output for user 
                    std::cout << " Feature " << counter_feature
                      << ": (" << final_feature.getRT()
                      << "," << final_feature.getMZ() << ") Qual.:"
                      << max_quality_feature << std::endl;
                    
                    // increase the total number of features
                    ++counter_feature;
       
					             // gather information for fitting summary
                    {
                      const Feature& f = this->features_->back();
  
					            // quality, correlation
                      CoordinateType corr = f.getOverallQuality();
                      summary.corr_mean += corr;
                      if (corr<summary.corr_min) summary.corr_min = corr;
                      if (corr>summary.corr_max) summary.corr_max = corr;

                      // charge
                      UInt ch = f.getCharge();
                      if (ch>= summary.charge.size())
                      {
                        summary.charge.resize(ch+1);
                      }
                      summary.charge[ch]++;

                      // MZ model type
                      const Param& p = f.getModelDescription().getParam();
                      ++summary.mz_model[ p.getValue("MZ") ];

                      // standard deviation of isotopic peaks
                      if (p.exists("MZ:isotope:stdev") && p.getValue("MZ:isotope:stdev")!=DataValue::EMPTY)
                      {
                        ++summary.mz_stdev[p.getValue("MZ:isotope:stdev")];
                      }
                    }
                    
                } // if
             
             }	// try
             catch(UnableToFit ex)
             {
								std::cout << "UnableToFit: " << ex.what() << std::endl;
 
                // set unused flag for all data points
                for (IndexSet::const_iterator it=region.begin(); it!=region.end(); ++it)
                {
                	this->ff_->getPeakFlag(*it) = UNUSED;
                }
              
                // gather information for fitting summary
                {
              	  ++summary.no_exceptions;
                  ++summary.exception[ex.getName()];
                }
    					} // catch
    	     	} // for (different charges ;-)
					} // for (boxes)
		
					this->ff_->endProgress();
              
         //---------------------------------------------------------------------------
         // print fitting summary
         //---------------------------------------------------------------------------
         {
            UInt size = this->features_->size();
            std::cout << size << " features were found. " << std::endl;

            // compute corr_mean
            summary.corr_mean /= size;

            std::cout << "FeatureFinder summary:\n"
                << "Correlation:\n\tminimum: " << summary.corr_min << "\n\tmean: " << summary.corr_mean
                << "\n\tmaximum: " << summary.corr_max << std::endl;

            std::cout << "Exceptions:\n";
            for (std::map<String,UInt>::const_iterator it=summary.exception.begin(); it!=summary.exception.end(); ++it)
            {
              std::cout << "\t" << it->first << ": " << it->second*100/summary.no_exceptions << "% (" << it->second << ")\n";
            }

            std::cout << "Chosen mz models:\n";
            for (std::map<String,UInt>::const_iterator it=summary.mz_model.begin(); it!=summary.mz_model.end(); ++it)
            {
              std::cout << "\t" << it->first << ": " << it->second*100/size << "% (" << it->second << ")\n";
            }

            std::cout << "Chosen mz stdevs:\n";
            for (std::map<float,UInt>::const_iterator it=summary.mz_stdev.begin(); it!=summary.mz_stdev.end(); ++it)
            {
              std::cout << "\t" << it->first << ": " << it->second*100/(size-summary.charge[0]) << "% (" << it->second << ")\n";
            }

            std::cout << "Charges:\n";
            for (UInt i=1; i<summary.charge.size(); ++i)
            {
              if (summary.charge[i]!=0)
              {
                std::cout << "\t+" << i << ": " << summary.charge[i]*100/(size-summary.charge[0]) << "% (" << summary.charge[i] << ")\n";
              }
            }
          }
            
          return;
          
        } // run

        static FeatureFinderAlgorithm<PeakType,FeatureType>* create()
        {
          return new FeatureFinderAlgorithmWavelet();
        }

        static const String getProductName()
        {
          return "isotope_wavelet";
        }
          
      protected:
        
        ///<Key: RT index, value: BoxElement_
        typedef typename IsotopeWaveletTransform<PeakType>::Box Box;
       
        /// The maximal charge state we will consider
        UInt max_charge_;
        /// The only parameter of the isotope wavelet
        CoordinateType ampl_cutoff_; 
        /// The number of susequent scans a pattern must cover in order to be considered as signal 
        UInt RT_votes_cutoff_; 
        /// The numer of scans we allow to be missed within RT_votes_cutoff_
        UInt RT_interleave_; 
        /// Negative or positive charged 
        Int mode_; 
        /// Charge threshold (in percent)
        CoordinateType charge_threshold_;

        virtual void updateMembers_()
        {
          max_charge_ = this->param_.getValue ("max_charge"); 
          ampl_cutoff_ = this->param_.getValue ("intensity_threshold");
          RT_votes_cutoff_ = this->param_.getValue ("rt_votes_cutoff");
          RT_interleave_ = this->param_.getValue ("rt_interleave");
          mode_ = this->param_.getValue ("recording_mode");
          IsotopeWavelet::setMaxCharge(max_charge_);
          charge_threshold_ = this->param_.getValue ("charge_threshold"); 
        }  
      
      private:
          
        /// Not implemented
        FeatureFinderAlgorithmWavelet& operator=(const FeatureFinderAlgorithmWavelet&);
        /// Not implemented
        FeatureFinderAlgorithmWavelet(const FeatureFinderAlgorithmWavelet&);

    };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMWAVELT_H
