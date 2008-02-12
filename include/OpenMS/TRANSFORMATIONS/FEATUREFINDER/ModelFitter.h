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


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MODELFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MODELFITTER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaIsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>
#include <OpenMS/MATH/STATISTICS/AsymmetricStatistics.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/CONCEPT/Factory.h>

#include <iostream>
#include <fstream>
#include <numeric>
#include <math.h>
#include <vector>
#include <set>

namespace OpenMS
{
    
    /**
        @brief Model fitter using gaussian or isotope model in mz and bigauss, lmagauss (bigauss with Levenberg-Marquardt aproximized parameters) or emg (exponent. modified Gaussian with lma aproximized arameters) in rt. For the isotope model different charges and deviations are tested.
    
        @ingroup FeatureFinder
    */
    template <class PeakType, class FeatureType> class ModelFitter :
    public FeaFiModule<PeakType,FeatureType>,
    public FeatureFinderDefs
    {
        public:
          
            /// IndexSet iterator
            typedef IndexSet::const_iterator IndexSetIter;
            /// Quality of a feature
            typedef Feature::QualityType QualityType;
            /// Single coordinate
            typedef Feature::CoordinateType CoordinateType;
            /// Single intensity 
            typedef Feature::IntensityType IntensityType;
            /// Isotope charge
            typedef Feature::ChargeType ChargeType;
            /// FeaFiModule
            typedef FeaFiModule<PeakType,FeatureType> Base;
            /// Raw data point type
            typedef RawDataPoint1D RawDataPointType;
            /// Raw data container type using for the temporary storage of the input data
            typedef DPeakArray<RawDataPointType > RawDataArrayType;
    
            enum 
            {
                RT = RawDataPoint2D::RT,
                MZ = RawDataPoint2D::MZ
            };

            /// Constructor
            ModelFitter(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff) :
                    Base(map,features,ff),
                    model2D_(),
                    mz_stat_(),
                    rt_stat_(),
                    monoisotopic_mz_( 0 ),
                    counter_( 1 ),
                    iso_stdev_first_( 0 ),
                    iso_stdev_last_( 0 ),
                    iso_stdev_stepsize_( 0 ),
                    first_mz_model_( 0 ),
                    last_mz_model_( 0 )
            {
                this->setName("ModelFitter");
                
                this->defaults_.setValue("fit_algorithm", "simple", "Fitting algorithm type (simplest, simple, advanced_isotope, ...).", false);
                
                this->defaults_.setValue( "max_iteration", 500, "Maximum number of iterations for fitting with Levenberg-marquardt algorithm.", true );
                this->defaults_.setValue( "deltaAbsError", 0.0001, "Absolute error used by the Levenberg-Marquardt algorithms.", true );
                this->defaults_.setValue( "deltaRelError", 0.0001, "Relative error used by the Levenberg-Marquardt algorithms.", true );
                
                this->defaults_.setValue( "tolerance_stdev_bounding_box", 3.0f, "Bounding box has range [minimim of data, maximum of data] enlarged by tolerance_stdev_bounding_box times the standard deviation of the data", true );
                
                this->defaults_.setValue( "intensity_cutoff_factor", 0.05f, "Cutoff peaks with a predicted intensity below intensity_cutoff_factor times the maximal intensity of the model", false );
                
                this->defaults_.setValue( "feature_intensity_sum", 1, "Determines what is reported as feature intensity.\n1: the sum of peak intensities;\n0: the maximum intensity of all peaks" , true);
                this->defaults_.setValue( "min_num_peaks:final", 5, "Minimum number of peaks left after cutoff. If smaller, feature will be discarded." , false);
                this->defaults_.setValue( "min_num_peaks:extended", 10, "Minimum number of peaks after extension. If smaller, feature will be discarded." , false);
                this->defaults_.setSectionDescription( "min_num_peaks", "Required number of peaks for a feature." );
                
                this->defaults_.setValue( "rt:interpolation_step", 0.2f, "Step size in seconds used to interpolate model for RT." , false);
                this->defaults_.setSectionDescription( "rt", "Model settings in RT dimension." );
                
                this->defaults_.setValue( "mz:interpolation_step", 0.03f, "Interpolation step size for m/z.", false );
                this->defaults_.setValue( "mz:model_type:first", 0, "Numeric id of first m/z model fitted (usually indicating the charge state), 0 = no isotope pattern (fit a single gaussian).", false );
                this->defaults_.setValue( "mz:model_type:last", 4, "Numeric id of last m/z model fitted (usually indicating the charge state), 0 = no isotope pattern (fit a single gaussian).", false );
                this->defaults_.setSectionDescription( "mz", "Model settings in m/z dimension." );
                
                this->defaults_.setValue( "quality:type", "Correlation_ModelFitter", "Type of the quality measure used to assess the fit of model vs data ('Correlation', 'EuclidianDistance', 'RankCorrelation').", true );
                this->defaults_.setValue( "quality:minimum", 0.65f, "Minimum quality of fit, features below this threshold are discarded." , false);
                this->defaults_.setSectionDescription( "quality", "Fitting quality settings." );
                
                this->defaults_.setValue( "isotope_model:stdev:first", 0.04f, "First standard deviation to be considered for isotope model.", true );
                this->defaults_.setValue( "isotope_model:stdev:last", 0.12f, "Last standard deviation to be considered for isotope model.", true );
                this->defaults_.setValue( "isotope_model:stdev:step", 0.04f, "Step size for standard deviations considered for isotope model.", true );
                this->defaults_.setSectionDescription( "isotope_model:stdev", "Instrument resolution settings for m/z dimension." );
                this->defaults_.setValue( "isotope_model:averagines:C", 0.0443f, "Number of C atoms per Dalton of the mass.", true );
                this->defaults_.setValue( "isotope_model:averagines:H", 0.007f, "Number of H atoms per Dalton of the mass.", true );
                this->defaults_.setValue( "isotope_model:averagines:N", 0.0012f, "Number of N atoms per Dalton of the mass.", true );
                this->defaults_.setValue( "isotope_model:averagines:O", 0.013f, "Number of O atoms per Dalton of the mass.", true );
                this->defaults_.setValue( "isotope_model:averagines:S", 0.00037f, "Number of S atoms per Dalton of the mass.", true);
                this->defaults_.setSectionDescription( "isotope_model:averagines", "Averagines are used to approximate the number of atoms (C,H,N,O,S) which a peptide of a given mass contains." ); 
                
                this->defaults_.setValue( "isotope_model:isotope:trim_right_cutoff", 0.001f, "Cutoff for averagine distribution, trailing isotopes below this relative intensity are not considered.", true );
                this->defaults_.setValue( "isotope_model:isotope:maximum", 100, "Maximum number of isotopes being used for the IsotopeModel.", true );
                this->defaults_.setValue( "isotope_model:isotope:distance", 1.000495f, "Distance between consecutive isotopic peaks.", true );
                this->defaults_.setSectionDescription( "isotope_model", "Settings of the isotope model (m/z)." );
    
                this->defaultsToParam_();
            }

            /// Destructor
            virtual ~ModelFitter()
            {
            }
            
            void setMonoIsotopicMass(CoordinateType mz)
            {
              monoisotopic_mz_ = mz;
            }
            
            /// Return next feature
            Feature fit(const ChargedIndexSet& index_set) throw (UnableToFit)
            {
                // Test the number of peaks (not enough peaks to fit)
                if ( index_set.size() < ( UInt ) ( this->param_.getValue( "min_num_peaks:extended" ) ) )
                {
                        String mess = String( "Skipping feature, IndexSet size too small: " ) + index_set.size();
                        throw UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-IndexSet", mess.c_str() );
                }
                
                // Calculate statistics for mz and rt
                mz_stat_.update
                (
                  Internal::IntensityIterator<ModelFitter>(index_set.begin(), this),
                  Internal::IntensityIterator<ModelFitter>(index_set.end(), this),
                  Internal::MzIterator<ModelFitter>(index_set.begin(), this)
                );
                rt_stat_.update
                (
                  Internal::IntensityIterator<ModelFitter>(index_set.begin(), this),
                  Internal::IntensityIterator<ModelFitter>(index_set.end(), this), 
                  Internal::RtIterator<ModelFitter>( index_set.begin(), this)
                );
                
                // Check charge estimate if charge is not specified by user
                if (index_set.charge_ != 0)
                {
                  first_mz_model_ = index_set.charge_;
                  last_mz_model_ = index_set.charge_;
                }
                
                std::cout << "Checking charge state from " << first_mz_model_ << " to " << last_mz_model_ << std::endl;
                
                // Compute model with the best correlation
                ProductModel<2>* final = 0;
                QualityType max_quality = fitLoop_(index_set, first_mz_model_, last_mz_model_, final);
                
                // model_desc.createModel() returns 0 if class model_desc is not initialized
                // in this case something went wrong during the model fitting and we stop.
                if ( ! final )
                {
                    throw UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-BadQuality", "Zero quality after fitting. Skipping this feature" );
                    delete final;
                }
                
                // find peak with highest predicted intensity to use as cutoff
                IntensityType model_max = 0;
                for ( IndexSetIter it = index_set.begin(); it != index_set.end(); ++it )
                {
                    IntensityType model_int = final->getIntensity( this->getPeakPos( *it ) );
                    if ( model_int > model_max ) model_max = model_int;
                }
                final->setCutOff( model_max * Real( this->param_.getValue( "intensity_cutoff_factor" ) ) );

                // Cutoff low intensities wrt to model maximum -> cutoff independent of scaling
                IndexSet model_set;
                for ( IndexSetIter it = index_set.begin(); it != index_set.end(); ++it )
                {
                    if ( final->isContained( this->getPeakPos( *it ) ) )
                    {
                        model_set.insert( *it );
                    }
                    else		// free dismissed peak via setting the appropriate flag
                    {
                        this->ff_->getPeakFlag( *it ) = UNUSED;
                    }
                }
                    
                // Print number of selected peaks after cutoff
                std::cout << " Selected " << model_set.size() << " from " << index_set.size() << " peaks.\n";

                // not enough peaks left for feature
                if ( model_set.size() < ( UInt ) ( this->param_.getValue( "min_num_peaks:final" ) ) )
                {
                    delete final;
                    throw UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__,"UnableToFit-FinalSet",String( "Skipping feature, IndexSet size after cutoff too small: " ) + model_set.size() );
                }

                std::vector<Real> data;
                data.reserve(model_set.size());
                std::vector<Real> model;
                model.reserve(model_set.size());

                for (IndexSet::iterator it=model_set.begin();it!=model_set.end();++it)
                {
                    data.push_back(this->getPeakIntensity(*it));
                    model.push_back(final->getIntensity(DPosition<2>(this->getPeakRt(*it),this->getPeakMz(*it))));
                }

                // fit has too low quality or fit was not possible i.e. because of zero stdev
                if ( max_quality < ( Real ) ( this->param_.getValue( "quality:minimum" ) ) )
                {
                    delete final;
                    String mess = String( "Skipping feature, correlation too small: " ) + max_quality;
                    throw UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-Correlation", mess.c_str() );
                }

                // Calculate intensity scaling
                IntensityType model_sum = 0;
                IntensityType data_sum = 0;
                IntensityType data_max = 0;
                for ( IndexSetIter it = model_set.begin(); it != model_set.end(); ++it )
                {
                    IntensityType model_int = final->getIntensity( this->getPeakPos( *it ) );
                    model_sum += model_int;
                    data_sum += this->getPeakIntensity( *it );
                    if ( this->getPeakIntensity( *it ) > data_max ) data_max = this->getPeakIntensity( *it );
                }
                
                // fit has too low quality or fit was not possible i.e. because of zero stdev
                if ( model_sum == 0 )
                {
                    delete final;
                    throw UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__,"UnableToFit-ZeroSum", "Skipping feature, model_sum zero." );
                }

                final->setScale( data_max / model_max );	// use max quotient instead of sum quotient

                // Build Feature
                // The feature coordinate in rt dimension is given
                // by the centroid of the rt model whereas the coordinate
                // in mz dimension is equal to the monoisotopic peak.
                Feature f;
                f.setModelDescription( ModelDescription<2>( final ) );
                f.setOverallQuality( max_quality );
                f.setRT( static_cast<InterpolationModel*>( final->getModel( RT ) ) ->getCenter() );
                f.setMZ( static_cast<InterpolationModel*>( final->getModel( MZ ) ) ->getCenter() );
                                                
                // feature charge ...	
                // if we used a simple Gaussian model to fit the feature, we can't say anything about
                // its charge state. The value 0 indicates that charge state is undetermined.
                if ( final->getModel( MZ ) ->getName() == "LmaIsotopeModel" )
                {
                    f.setCharge( static_cast<LmaIsotopeModel*>( final->getModel( MZ ) ) ->getCharge() );
                }
                else if (final->getModel( MZ ) ->getName() == "IsotopeModel")
                {
                    f.setCharge( static_cast<IsotopeModel*>( final->getModel( MZ ) ) ->getCharge() );
                }
                else
                {
                    f.setCharge( 0 );
                }
            
                // feature intensity
                Int const intensity_choice = this->param_.getValue( "feature_intensity_sum" );
                IntensityType feature_intensity = 0.0;
                if ( intensity_choice == 1 )
                {
                    // intensity of the feature is the sum of all included data points
                    for ( IndexSetIter it = model_set.begin(); it != model_set.end(); ++it )
                    {
                        feature_intensity += this->getPeakIntensity( *it );
                    }
                }
                else
                {
                    // feature intensity is the maximum intensity of all peaks
                    for ( IndexSetIter it = model_set.begin(); it != model_set.end(); ++it )
                    {
                        if ( this->getPeakIntensity( *it ) > feature_intensity )
                        {
                                feature_intensity = this->getPeakIntensity( *it );
                        }
                    }
                }
                              
                f.setIntensity( feature_intensity );
                this->addConvexHull( model_set, f );

                std::cout << __FILE__ << ':' << __LINE__ << ": " << QDateTime::currentDateTime().toString( "yyyy-MM-dd hh:mm:ss" ).toStdString() << " Feature " << counter_
                    << ": (" << f.getRT()
                    << "," << f.getMZ() << ") Qual.:"
                    << max_quality << std::endl;

                // Calculate quality	
                //RT fit
                data.clear();
                model.clear();
                for (IndexSet::iterator it=model_set.begin();it!=model_set.end();++it)
                {
                    data.push_back(this->getPeakIntensity(*it));
                    model.push_back((final->getModel(RT))->getIntensity(this->getPeakRt(*it)));
                }
                f.setQuality( RT, Math::BasicStatistics<Real>::pearsonCorrelationCoefficient(data.begin(), data.end(), model.begin(), model.end()));
                
                //MZ fit
                data.clear();
                model.clear();
                for (IndexSet::iterator it=model_set.begin();it!=model_set.end();++it)
                {
                    data.push_back(this->getPeakIntensity(*it));
                    model.push_back((final->getModel(MZ))->getIntensity(this->getPeakMz(*it)));
                }
                f.setQuality( MZ, Math::BasicStatistics<Real>::pearsonCorrelationCoefficient(data.begin(), data.end(), model.begin(), model.end()));

                // Save meta data in feature for TOPPView
                std::stringstream meta ;
                meta << "Feature #" << counter_ << ", +"	<< f.getCharge() << ", " << index_set.size() << "->" << model_set.size()
                        << ", Corr: (" << max_quality << "," << f.getQuality( RT ) << "," << f.getQuality( MZ ) << ")";
                f.setMetaValue( 3, String( meta.str() ) );

                // Count features
                ++counter_;
                
                delete final;

                return f;
            }
            
    protected:

            virtual void updateMembers_()
            {
                algorithm_ = this->param_.getValue( "fit_algorithm" );
                
                max_iteration_ = this->param_.getValue("max_iteration");
                deltaAbsError_ = this->param_.getValue("deltaAbsError");
                deltaRelError_ = this->param_.getValue("deltaRelError");
                
                tolerance_stdev_box_ = this->param_.getValue( "tolerance_stdev_bounding_box" );
                
                interpolation_step_mz_ = this->param_.getValue( "mz:interpolation_step" );
                interpolation_step_rt_ = this->param_.getValue( "rt:interpolation_step" );
    
                iso_stdev_first_ = this->param_.getValue( "isotope_model:stdev:first" );
                iso_stdev_last_ = this->param_.getValue( "isotope_model:stdev:last" );
                iso_stdev_stepsize_ = this->param_.getValue( "isotope_model:stdev:step" );
    
                first_mz_model_ = ( Int ) this->param_.getValue( "mz:model_type:first" );
                last_mz_model_ = ( Int ) this->param_.getValue( "mz:model_type:last" );
            }
             
           /// main fit loop
           QualityType fitLoop_(const ChargedIndexSet& set, Int& first_mz, Int& last_mz, ProductModel<2>*& final) 
           {
             // Projection    
              doProjectionDim_(set, rt_input_data_, RT, algorithm_);
              total_intensity_mz_ = doProjectionDim_(set, mz_input_data_, MZ, algorithm_);
            
              // Fit rt model
              fitDim_(RT, algorithm_);
              
              // Fit mz model ... test different charge states and stdevs
              QualityType quality_mz = 0.0;
              QualityType max_quality_mz = -std::numeric_limits<QualityType>::max();
      
              std::map<QualityType,ProductModel<2> > model_map;
              for ( Real stdev = iso_stdev_first_; stdev <= iso_stdev_last_; stdev += iso_stdev_stepsize_)
              {
                  for (Int mz_fit_type = first_mz; mz_fit_type <= last_mz; ++mz_fit_type)
                  {
                    charge_ = mz_fit_type;
                    isotope_stdev_ = stdev;
                    quality_mz =  fitDim_(MZ, algorithm_);
                                      
                    if (quality_mz > max_quality_mz)
                    {
                        max_quality_mz = quality_mz;
                        model_map.insert( std::make_pair( quality_mz, model2D_) ); 
                    }
                  }
              }
              
              std::map<QualityType,ProductModel<2> >::iterator it_map = model_map.find(max_quality_mz);
              final = new ProductModel<2>((*it_map).second);
              
              // Compute overall quality
              QualityType max_quality = 0.0;
              max_quality = evaluate_(set, final, algorithm_);
                                        
              return max_quality;
            }
            
            /// evaluate 2d-model
            QualityType evaluate_(const IndexSet& set, ProductModel<2>*& final, String algorithm)
            {
                QualityType quality = 1.0;
              
                // Calculate the pearson correlation coefficient for the values in [begin_a, end_a) and [begin_b, end_b)
                if (algorithm!="")
                {
                    std::vector<Real> real_data;
                    real_data.reserve(set.size());
                    std::vector<Real> model_data;
                    model_data.reserve(set.size());
          
                    for (IndexSet::iterator it=set.begin(); it != set.end(); ++it)
                    {
                      real_data.push_back(this->getPeakIntensity(*it));
                      model_data.push_back(final->getIntensity(DPosition<2>(this->getPeakRt(*it),this->getPeakMz(*it))));
                    }
                    
                    quality = basic_stat_.pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());
                }
              
                if (isnan(quality)) quality = -1.0; 
    
                return quality;
            }
            
            /// do a fit for a given dimension and algorithm
            QualityType fitDim_(Int dim, String algorithm)  
            {
              QualityType quality;
              Param param;
              Fitter1D* fitter;
              InterpolationModel* model = 0;
              
              if (dim==RT)
              {
                if (algorithm=="simplest")
                {
                    // Parameter for rt (BiGauss)
                    param.setValue( "tolerance_stdev_bounding_box", tolerance_stdev_box_);
                    param.setValue( "statistics:mean", rt_stat_.mean() );
                    param.setValue( "statistics:variance", rt_stat_.variance() );
                    param.setValue( "statistics:variance1", rt_stat_.variance1() );
                    param.setValue( "statistics:variance2", rt_stat_.variance2() );
                    param.setValue( "interpolation_step", interpolation_step_rt_ );
                    
                    fitter = Factory<Fitter1D >::create("BiGaussFitter1D");
                }
                else
                {
                    // Parameter for rt (EMG)
                    param.setValue( "tolerance_stdev_bounding_box", tolerance_stdev_box_);
                    param.setValue( "statistics:mean", rt_stat_.mean() );
                    param.setValue( "statistics:variance", rt_stat_.variance() );
                    param.setValue( "interpolation_step", interpolation_step_rt_ );
                    param.setValue( "max_iteration", max_iteration_);
                    param.setValue( "deltaAbsError", deltaAbsError_);
                    param.setValue( "deltaRelError", deltaRelError_);
                      
                    fitter = Factory<Fitter1D >::create("EmgFitter1D");
                }

                // Set parameter for fitter                
                fitter->setParameters( param );
      
                // Construct model for rt
                quality = fitter->fit1d(rt_input_data_, model);
              }
              else
              {
                  // Parameter for mz (Isotope)
                  param.setValue( "tolerance_stdev_bounding_box", tolerance_stdev_box_);
                  param.setValue( "statistics:mean", mz_stat_.mean() );
                  param.setValue( "statistics:variance", mz_stat_.variance() );
                  param.setValue( "interpolation_step", interpolation_step_mz_ );
                  param.setValue( "charge", charge_ );
                  param.setValue( "isotope:stdev", isotope_stdev_ );
                  
                  // monoisotopic mass is known :-)
                  if (monoisotopic_mz_ != 0)
                  {
                    param.setValue( "statistics:mean", monoisotopic_mz_ );
                  }
                  
                  if (algorithm=="lmaiso")
                  {
                    param.setValue( "total_intensity", total_intensity_mz_ );
                    param.setValue( "max_iteration", max_iteration_);
                    param.setValue( "deltaAbsError", deltaAbsError_);
                    param.setValue( "deltaRelError", deltaRelError_);
                    
                     // monoisotopic mass is known :-)
                    if (monoisotopic_mz_ != 0)
                    {
                      param.setValue( "monoisotopic_mass", monoisotopic_mz_ );
                    }
                 
                    fitter = Factory<Fitter1D >::create("LmaIsotopeFitter1D");
                  }
                  else
                  {
                    fitter = Factory<Fitter1D >::create("IsotopeFitter1D");
                  }
                  
                  fitter->setParameters( param );
                  
                  // Construct model for mz
                  quality = fitter->fit1d(mz_input_data_, model);
              }

              // Check quality
              if (isnan(quality) ) quality = -1.0;
              
              // Set model in 2D-model						
              model2D_.setModel(dim, model);
              
              delete(fitter);
              
              return quality;
            }
           
            /// copy the raw data into 1-dim. DPeakArray 
            CoordinateType doProjectionDim_(const ChargedIndexSet& index_set, RawDataArrayType& set, Int dim, String algorithm)
            {
              CoordinateType total_intensity = 0;
         
              if (algorithm!="")
              {
                std::map<CoordinateType,CoordinateType> data_map;
                 
                if (dim==MZ)
                {
                  for ( IndexSet::const_iterator it = index_set.begin(); it != index_set.end(); ++it)
                  {   
                    data_map[this->getPeakMz(*it)] += this->getPeakIntensity(*it);
                  }  
                }
                else
                {
                  for ( IndexSet::const_iterator it = index_set.begin(); it != index_set.end(); ++it)
                  { 
                    data_map[this->getPeakRt(*it)] += this->getPeakIntensity(*it);
                  }  
                }
                 
                // Copy the raw data into a DPeakArray<DRawDataPoint<D> >
                set.resize(data_map.size());
                std::map<CoordinateType,CoordinateType>::iterator it;
                UInt i=0;
                for ( it=data_map.begin() ; it != data_map.end(); ++it, ++i )
                {
                  set[i].setPosition((*it).first);
                  set[i].setIntensity((*it).second);
                  total_intensity += (*it).second;
                }   
                
                data_map.clear();
              }            
              
              return total_intensity;
            }
            
            /// 2D model
            ProductModel<2> model2D_;
            /// statistics for mz
            Math::BasicStatistics<> mz_stat_;
            /// statistics for rt
            Math::AsymmetricStatistics<> rt_stat_;
            /// mz raw data
            RawDataArrayType mz_input_data_;
            /// rt raw data
            RawDataArrayType rt_input_data_;
            /// tolerance used for bounding box
            CoordinateType tolerance_stdev_box_;
            /// monoistopic mass
            CoordinateType monoisotopic_mz_;
            /// counts features (used for debug output only)
            UInt counter_;
            /// interpolation step size (in m/z)
            CoordinateType interpolation_step_mz_;
            /// interpolation step size (in retention time)
            CoordinateType interpolation_step_rt_;
            /// first stdev
            CoordinateType iso_stdev_first_;
            /// last stdev
            CoordinateType iso_stdev_last_;
            /// step size
            CoordinateType iso_stdev_stepsize_;
            /// first mz model (0=Gaussian, 1....n = charge )
            Int first_mz_model_;			
            /// last mz model
            Int last_mz_model_;
            /// isotope charge
            ChargeType charge_;
            /// isotope stdev
            CoordinateType isotope_stdev_;
            /// algorithm
            String algorithm_;
            /// Maximum number of iterations
            Int max_iteration_;
            /** Test for the convergence of the sequence by comparing the last iteration step dx with the absolute error epsabs and relative error epsrel to the current position x */
            /// Absolute error
            CoordinateType deltaAbsError_;
            /// Relative error
            CoordinateType deltaRelError_;
            /// statistics
            Math::BasicStatistics<> basic_stat_;
            /// area under mz curve
            CoordinateType total_intensity_mz_;

    private:

            /// Not implemented
            ModelFitter();
            /// Not implemented
            ModelFitter& operator=(const ModelFitter&);
            /// Not implemented
            ModelFitter(const ModelFitter&);

	};

}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_MODELFITTER_H
