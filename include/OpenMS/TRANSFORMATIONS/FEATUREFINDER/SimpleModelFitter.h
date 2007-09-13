// -*- Mode: C++; tab-width: 2; -*-
// 
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
// $Maintainer: Clemens Groepl, Marcel Grunert $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEMODELFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEMODELFITTER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/AsymmetricStatistics.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LogNormalModel.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/SYSTEM/StopWatch.h>

#include <iostream>
#include <fstream>
#include <numeric>
#include <math.h>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

namespace OpenMS
{

	namespace Internal
	{
		// Helper struct for SimpleModelFitter
		struct ExpFitPolyData
		{
			size_t n;
			std::string profile;
		};
	}

	/**
		@brief Extended model fitter using gaussian or isotope model in mz and
		bigauss, lmagauss (bigauss with Levenberg-Marquardt aproximized
		parameters) or emg (exponent. modified Gaussian with lma aproximized
		parameters) in rt.

		For the isotope model different charges and deviations are tested.

		@ref SimpleModelFitter_Parameters are explained on a separate page.

		@todo Check use of Enums for RT and m/z fit. They destroy the factory concept! (Clemens)

		@ingroup FeatureFinder
		*/
	template <class PeakType, class FeatureType>
		class SimpleModelFitter
		: public FeaFiModule<PeakType,FeatureType>,
		public FeatureFinderDefs
	{
		public:
			///
			typedef FeaFiModule<PeakType,FeatureType> Base;

			/// Ion count
			typedef typename Base::IntensityType IntensityType;
			/// Quality of a feature
			typedef Feature::QualityType QualityType;
			/// IndexSet iterator
			typedef IndexSet::const_iterator IndexSetIter;
			/// Single coordinate
			typedef typename Base::CoordinateType Coordinate;
			///	Single coordinate
			typedef Feature::CoordinateType CoordinateType;

			enum RtFitting{ RTGAUSS=0, LMAGAUSS=1, EMGAUSS=2, BIGAUSS=3, LOGNORMAL=4 };
			enum MzFitting{ MZGAUSS=0, CHARGE1=1, CHARGE2=2, CHARGE3=3, CHARGE4=4	};

			enum 
			{
				RT = RawDataPoint2D::RT,
				MZ = RawDataPoint2D::MZ
			};

			/// Constructor
			SimpleModelFitter(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff)
				: Base(map,features,ff),
				model2D_(),
				mz_stat_(),
				rt_stat_(),
				stdev_mz_( 0 ),
				stdev_rt1_( 0 ),
				stdev_rt2_( 0 ),
				min_(),
				max_(),
				counter_( 1 ),
				iso_stdev_first_( 0 ),
				iso_stdev_last_( 0 ),
				iso_stdev_stepsize_( 0 ),
				first_mz_model_( 0 ),
				last_mz_model_( 0 )
		{
			this->setName("SimpleModelFitter");

			this->defaults_.setValue( "tolerance_stdev_bounding_box", 3.0f, "Bounding box has range [minimim of data, maximum of data] enlarged by tolerance_stdev_bounding_box times the standard deviation of the data" );
			this->defaults_.setValue( "intensity_cutoff_factor", 0.05f, "Cutoff peaks with a predicted intensity below intensity_cutoff_factor times the maximal intensity of the model" );
			this->defaults_.setValue( "feature_intensity_sum", 1, "Determines what is reported as feature intensity.\n1: the sum of peak intensities;\n0: the maximum intensity of all peaks" );

			this->defaults_.setValue( "min_num_peaks:final", 5, "Minimum number of peaks left after cutoff. If smaller, feature will be discarded." );
			this->defaults_.setValue( "min_num_peaks:extended", 10, "Minimum number of peaks after extension. If smaller, feature will be discarded." );
			this->defaults_.setSectionDescription( "min_num_peaks", "Required number of peaks for a feature." );

			this->defaults_.setValue( "rt:interpolation_step", 0.2f, "Step size in seconds used to interpolate model for RT." );
			this->defaults_.setValue( "rt:max_iteration", 500, "Maximum number of iterations for RT fitting." );
			this->defaults_.setValue( "rt:deltaAbsError", 0.0001, "Absolute error used by the Levenberg-Marquardt algorithms." );
			this->defaults_.setValue( "rt:deltaRelError", 0.0001, "Relative error used by the Levenberg-Marquardt algorithms." );
			this->defaults_.setValue( "rt:profile", "EMG", "Type of RT model. Possible models are 'LmaGauss', 'EMG' and 'LogNormal'." );
			this->defaults_.setSectionDescription( "rt", "Model settings in RT dimension." );

			this->defaults_.setValue( "mz:interpolation_step", 0.03f, "Interpolation step size for m/z." );
			this->defaults_.setValue( "mz:model_type:first", 0, "Numeric id of first m/z model fitted (usually indicating the charge state), 0 = no isotope pattern (fit a single gaussian)." );
			this->defaults_.setValue( "mz:model_type:last", 4, "Numeric id of last m/z model fitted (usually indicating the charge state), 0 = no isotope pattern (fit a single gaussian)." );
			this->defaults_.setSectionDescription( "mz", "Model settings in m/z dimension." );

			this->defaults_.setValue( "quality:type", "Correlation", "Type of the quality measure used to assess the fit of model vs data ('Correlation','EuclidianDistance','RankCorrelation')." );
			this->defaults_.setValue( "quality:minimum", 0.65f, "Minimum quality of fit, features below this threshold are discarded." );
			this->defaults_.setSectionDescription( "quality", "Fitting quality settings." );

			this->defaults_.setValue( "isotope_model:stdev:first", 0.04f, "First standard deviation to be considered for isotope model." );
			this->defaults_.setValue( "isotope_model:stdev:last", 0.12f, "Last standard deviation to be considered for isotope model." );
			this->defaults_.setValue( "isotope_model:stdev:step", 0.04f, "Step size for standard deviations considered for isotope model." );
			this->defaults_.setSectionDescription( "isotope_model:stdev", "Instrument resolution settings for m/z dimension." );

			this->defaults_.setValue( "isotope_model:averagines:C", 0.0443f, "Number of C atoms per Dalton of the mass." );
			this->defaults_.setValue( "isotope_model:averagines:H", 0.007f, "Number of H atoms per Dalton of the mass." );
			this->defaults_.setValue( "isotope_model:averagines:N", 0.0012f, "Number of N atoms per Dalton of the mass." );
			this->defaults_.setValue( "isotope_model:averagines:O", 0.013f, "Number of O atoms per Dalton of the mass." );
			this->defaults_.setValue( "isotope_model:averagines:S", 0.00037f, "Number of S atoms per Dalton of the mass." );
			this->defaults_.setSectionDescription( "isotope_model:averagines", "Averagines are used to approximate the number of atoms (C,H,N,O,S) which a peptide of a given mass contains." );

			this->defaults_.setValue( "isotope_model:isotope:trim_right_cutoff", 0.001f, "Cutoff for averagine distribution, trailing isotopes below this relative intensity are not considered." );
			this->defaults_.setValue( "isotope_model:isotope:maximum", 100, "Maximum number of isotopes being used for the IsotopeModel." );
			this->defaults_.setValue( "isotope_model:isotope:distance", 1.000495f, "Distance between consecutive isotopic peaks." );
			this->defaults_.setSectionDescription( "isotope_model", "Settings of the isotope model (m/z)." );

			this->defaultsToParam_();
		}

			/// Destructor
			virtual ~SimpleModelFitter()
			{
			}

			/// Return next feature
			Feature fit(const ChargedIndexSet& index_set) throw (UnableToFit)
			{
				// not enough peaks to fit
				if ( index_set.size() < ( UInt ) ( this->param_.getValue( "min_num_peaks:extended" ) ) )
				{
					String mess = String( "Skipping feature, IndexSet size too small: " ) + index_set.size();
					throw UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-IndexSet", mess.c_str() );
				}

				ModelDescription<2> model_desc;
				QualityType quality = 0.0;
				QualityType max_quality = -std::numeric_limits<double>::max();

				// Calculate statistics
				mz_stat_.update(Internal::IntensityIterator<SimpleModelFitter>(index_set.begin(), this), Internal::IntensityIterator<SimpleModelFitter>(index_set.end(), this), Internal::MzIterator<SimpleModelFitter>(index_set.begin(), this));
				rt_stat_.update(Internal::IntensityIterator<SimpleModelFitter>(index_set.begin(), this), Internal::IntensityIterator<SimpleModelFitter>(index_set.end(), this), Internal::RtIterator<SimpleModelFitter>( index_set.begin(), this));

				// Calculate bounding box
				IndexSetIter it = index_set.begin();
				min_ = max_ = this->getPeakPos( *it );
				for ( ++it; it != index_set.end(); ++it )
				{
					CoordinateType tmp = this->getPeakMz( *it );
					if ( min_[ MZ ] > tmp ) min_[ MZ ] = tmp;
					if ( max_[ MZ ] < tmp ) max_[ MZ ] = tmp;
					tmp = this->getPeakRt( *it );
					if ( min_[ RT ] > tmp ) min_[ RT ] = tmp;
					if ( max_[ RT ] < tmp ) max_[ RT ] = tmp;
				}

				// Enlarge the bounding box by a few multiples of the standard deviation
				{
					double const tolerance_stdev_box = this->param_.getValue( "tolerance_stdev_bounding_box" );
					stdev_mz_ = sqrt ( mz_stat_.variance() ) * tolerance_stdev_box;
					min_[ MZ ] -= stdev_mz_;
					max_[ MZ ] += stdev_mz_;

					stdev_rt1_ = sqrt ( rt_stat_.variance1() ) * tolerance_stdev_box;
					stdev_rt2_ = sqrt ( rt_stat_.variance2() ) * tolerance_stdev_box;
					min_[ RT ] -= stdev_rt1_;
					max_[ RT ] += stdev_rt2_;
				}

				// Create a vector with RT-values and intensity; compute the parameters
				// (intial values) for the EMG and Gauss function; and finally, optimize
				// the parameters with Levenberg-Marquardt algorithms
				if ( profile_ == "LmaGauss" || profile_ == "EMG" || profile_ == "LogNormal" )
				{
					setInitialParameters( index_set );
					if ( symmetric_ == false ) optimize();

					if ( gsl_status_ != "success" )
					{
						std::cout << profile_ + " status: " + gsl_status_ << std::endl;
						//throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,"UnableToFit-BadQuality",String("Skipping feature, " + profile_ + " status: " + gsl_status_));
					}
				}

				// IWASHERE;

				/// Test different charges and stdevs
				Int first_mz = first_mz_model_;
				Int last_mz = last_mz_model_;

				/// Check charge estimate if charge is not specified by user
				if ( index_set.charge_ != 0 /*&& (iso_stdev_first_ != iso_stdev_last_)*/ )
				{
					// 			first_mz = index_set.charge_;
					// 			last_mz = index_set.charge_;
					// 	first_mz = (index_set.charge_ - 1);
					//	last_mz = (index_set.charge_ + 1);
				}
				std::cout << "Checking charge state from " << first_mz << " to " << last_mz << std::endl;

				// IWASHERE;

				ProductModel<2>* final = 0;	// model  with best correlation

				for ( float stdev = iso_stdev_first_; stdev <= iso_stdev_last_; stdev += iso_stdev_stepsize_ )
				{
					for ( Int mz_fit_type = first_mz; mz_fit_type <= last_mz; ++mz_fit_type )
					{
						// IWASHEREMSG(profile_);
						if ( profile_ == "LmaGauss" )
						{
							quality = fit_( index_set, static_cast<MzFitting>( mz_fit_type ), LMAGAUSS, stdev );
						}
						else if ( profile_ == "EMG" && symmetric_ == false )
						{
							quality = fit_( index_set, static_cast<MzFitting>( mz_fit_type ), EMGAUSS, stdev );
						}
						else if ( profile_ == "LogNormal" && symmetric_ == false && symmetry_ != 1 && symmetry_ != 0 )
						{
							quality = fit_( index_set, static_cast<MzFitting>( mz_fit_type ), LOGNORMAL, stdev );
						}
						else
						{
							quality = fit_( index_set, static_cast<MzFitting>( mz_fit_type ), BIGAUSS, stdev );
						}
						
						// IWASHERE;

						if ( quality > max_quality )
						{
							max_quality = quality;
							//model_desc = ModelDescription<2>(&model2D_);
							final = new ProductModel<2>( model2D_ );	// store model
						}

					}
				// IWASHERE;
				}

				// IWASHERE;
				
				// model with highest correlation
				//ProductModel<2>* final = dynamic_cast< ProductModel<2>* >(model_desc.createModel());

				// model_desc.createModel() returns 0 if class model_desc is not initialized
				// in this case something went wrong during the modelfitting and we stop.
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
				final->setCutOff( model_max * float( this->param_.getValue( "intensity_cutoff_factor" ) ) );

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
					throw UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__,
							"UnableToFit-FinalSet",
							String( "Skipping feature, IndexSet size after cutoff too small: " ) + model_set.size() );
				}

				std::vector<Real> data(model_set.size());
				std::vector<Real> model(model_set.size());
				
				for (IndexSet::iterator it=model_set.begin();it!=model_set.end();++it)
				{
					data.push_back(this->getPeakIntensity(*it));
					model.push_back(final->getIntensity(DPosition<2>(this->getPeakRt(*it),this->getPeakMz(*it))));
				}

				max_quality = Math::BasicStatistics<Real>::pearsonCorrelationCoefficient(data.begin(), data.end(), model.begin(), model.end());

				// fit has too low quality or fit was not possible i.e. because of zero stdev
				if ( max_quality < ( float ) ( this->param_.getValue( "quality:minimum" ) ) )
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
					throw UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__,
							"UnableToFit-ZeroSum", "Skipping feature, model_sum zero." );
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
				if ( final->getModel( MZ ) ->getName() == "IsotopeModel" )
				{
					f.setCharge( static_cast<IsotopeModel*>( final->getModel( MZ ) ) ->getCharge() );
				}
				// if we used a simple Gaussian model to fit the feature, we can't say anything about
				// its charge state. The value 0 indicates that charge state is undetermined.
				else
				{
					f.setCharge( 0 );
				}

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
							feature_intensity = this->getPeakIntensity( *it );
					}
				}

				f.setIntensity( feature_intensity );
				this->addConvexHull( model_set, f );

				std::cout << QDateTime::currentDateTime().toString( "yyyy-MM-dd hh:mm:ss" ).toStdString() << " Feature " << counter_
					<< ": (" << f.getRT()
					<< "," << f.getMZ() << ") Qual.:"
					<< max_quality << "\n";

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

				// save meta data in feature for TOPPView
				std::stringstream meta ;
				meta << "Feature #" << counter_ << ", +"	<< f.getCharge() << ", " << index_set.size() << "->" << model_set.size()
					<< ", Corr: (" << max_quality << "," << f.getQuality( RT ) << "," << f.getQuality( MZ ) << ")";
				f.setMetaValue( 3, String( meta.str() ) );


#ifdef DEBUG_FEATUREFINDER
				std::cout << "Feature charge: " << f.getCharge() << std::endl;
				std::cout << "Feature quality in mz: " << f.getQuality( MZ ) << std::endl;
#endif

#ifdef DEBUG_FEATUREFINDER
				// write debug output
				CoordinateType rt = f.getRT();
				CoordinateType mz = f.getMZ();

				// write feature model
				String fname = String( "model" ) + counter_ + "_" + rt + "_" + mz;
				ofstream file( fname.c_str() );
				for ( IndexSetIter it = model_set.begin(); it != model_set.end(); ++it )
				{
					DPosition<2> pos = this->getPeakPos( *it );
					if ( final->isContained( pos ) )
					{
						file << pos[ RT ] << " " << pos[ MZ ] << " " << final->getIntensity( this->getPeakPos( *it ) ) << "\n";
					}
				}
				file.close();

				// wrote peaks remaining after model fit
				fname = String( "feature" ) + counter_ + "_" + rt + "_" + mz;
				ofstream file2( fname.c_str() );
				for ( IndexSetIter it = model_set.begin(); it != model_set.end(); ++it )
				{
					DPosition<2> pos = this->getPeakPos( *it );
					if ( final->isContained( pos ) )
					{
						file2 << pos[ RT ] << " " << pos[ MZ ] << " " << this->getPeakIntensity( *it ) << "\n";
					}
				}
				file2.close();
#endif

				++counter_;

				delete final;

				return f;
			} // fit()
			

			/// create a vector with RT-values & Intensities and compute the parameters (intial values) for the EMG, Gauss and logNormal function
			void setInitialParameters (const IndexSet& set)
			{
				// sum over all intensities
				double sum = 0.0;

				// iterate over all points of the signal
				for ( IndexSet::const_iterator it = set.begin(); it != set.end(); it++ )
				{
					// store the current rt-position and signal
					float position = this->getPeakRt( *it );
					float signal = this->getPeakIntensity( *it );

					//float mz = this->getPeakMz(*it);
					sum += signal;

					//orgFile << position << "  " << mz << " " << signal << "\n";

					// fill vectors with rt-postion and signal
					if ( positionsDC_.empty() || positionsDC_.back() != position )
					{
						positionsDC_.push_back( position );
						signalDC_.push_back( signal );
					}
					else
					{
						signal += signalDC_.back();
						signalDC_.pop_back();
						signalDC_.push_back( signal );
					}
				}

				// calculate the median
				int median = 0;
				float count = 0.0;
				for ( size_t current_point = 0; current_point < positionsDC_.size();current_point++ )
				{
					count += signalDC_[ current_point ];
					if ( count <= sum / 2 )
					{
						median = current_point;
					}
				}

				double sumS = 0.0;
				for ( size_t current_point = 0; current_point < positionsDC_.size();current_point++ )
				{
					sumS += pow( ( positionsDC_[ current_point ] - positionsDC_[ median ] ), 2 );
				}

				// calculate the stardard deviation
				standard_deviation_ = sqrt( sumS / ( positionsDC_.size() - 1 ) );

				// set expeceted value
				expected_value_ = positionsDC_[ median ]; //rt_stat_.mean();

				// calculate the heigth of the peak
				height_ = signalDC_[ median ];

				// calculate the width of the peak
				// rt-values with intensity zero are not allowed for calculation of the width
				width_ = fabs( positionsDC_[ positionsDC_.size() - 1 ] - positionsDC_[ 0 ] );

				// calculate retention time
				retention_ = positionsDC_[ median ];

				// default is an asymmetric peak
				symmetric_ = false;

				// calculate the symmetry (fronted peak: s<1 , tailed peak: s>1)
				symmetry_ = fabs( positionsDC_.back() - positionsDC_[ median ] ) / fabs( positionsDC_[ median ] - positionsDC_.front() );

				// check the symmetry
				if ( isinf( symmetry_ ) || isnan( symmetry_ ) )
				{
					symmetric_ = true;
					symmetry_ = 10;
				}

				// optimize the symmetry
				if ( profile_ == "LogNormal" )
				{
					// The computations can lead to an overflow error at very low values of symmetry (s~0).
					if ( symmetry_ <= 0.8 ) symmetry_ = 0.8;
					if ( symmetry_ == 1 ) symmetry_ = 1.1;
					if ( symmetry_ >= 1.5 ) symmetry_ = 1.4;

					// it is better to proceed from narrow peaks
					width_ /= 2;
				}
				else
				{
					// The computations can lead to an overflow error at very low values of symmetry (s~0). For s~5 the parameter can be aproximized by the Levenberg-Marquardt argorithms. (the other parameters are much greater than one)
					if ( symmetry_ < 1 ) symmetry_ += 5;

					// it is better for the emg function to proceed from narrow peaks
					width_ = symmetry_;
				}

				/* set the parameter r of the log normal function;
					 r is the ratio between h and the height at which w and s are computed;
					 r = 2, see "Mathematical functions for representation of chromatographic peaks", V.B. Di Marco(2001)
					 */
				r_ = 2;
			}

			/// perform a nonlinear optimization
			void optimize()
			{
				const gsl_multifit_fdfsolver_type * T;
				gsl_multifit_fdfsolver *s;

				int status;
				size_t iter = 0;
				const size_t n = positionsDC_.size();

				// number of parameter to be optimize
				unsigned int p = 0;
				if ( profile_ == "LmaGauss" ) p = 3;
				else if ( profile_ == "LogNormal" ) p = 4; //5;
				else p = 4;

				// gsl always excepts N>=p or default gsl error handler invoked, cause Jacobian be rectangular M x N with M>=N
				if ( n < p ) throw UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-FinalSet", "Skipping feature, gsl always expects N>=p" );

				gsl_matrix *covar = gsl_matrix_alloc( p, p );
				gsl_multifit_function_fdf f;

				double x_init_normal[ 3 ] = { standard_deviation_, expected_value_, height_ };
				double x_init_emg[ 4 ] = { height_, width_, symmetry_, retention_ };
				//double x_init_lognormal[5] = { height_, width_, symmetry_, retention_, r_ };
				double x_init_lognormal[ 4 ] = { height_, width_, symmetry_, retention_ };

				gsl_vector_view x;

				if ( profile_ == "LmaGauss" )
				{
					x = gsl_vector_view_array( x_init_normal, p );
				}
				else
				{
					if ( profile_ == "LogNormal" )
					{
						x = gsl_vector_view_array( x_init_lognormal, p );
					}
					else
					{
						x = gsl_vector_view_array( x_init_emg, p );
					}
				}

				const gsl_rng_type * type;
				gsl_rng * r;
				gsl_rng_env_setup();
				type = gsl_rng_default;
				r = gsl_rng_alloc ( type );

				struct Internal::ExpFitPolyData d =
				{
					n, profile_
				};
				///@todo Fix this. I have no clue how GSL works... (Clemens,Marcel)
				f.f = &(SimpleModelFitter::residualDC);
				f.df = &(SimpleModelFitter::jacobianDC);
				f.fdf = &(SimpleModelFitter::evaluateDC);
				f.n = n;
				f.p = p;
				f.params = &d;

				T = gsl_multifit_fdfsolver_lmsder;
				s = gsl_multifit_fdfsolver_alloc( T, n, p );
				gsl_multifit_fdfsolver_set( s, &f, &x.vector );

#ifdef DEBUG_FEATUREFINDER
				if ( profile_ == "LmaGauss" )
				{
					printf ( "before loop iter: %4u x = % 15.8f % 15.8f % 15.8f |f(x)| = %g\n", iter,
							gsl_vector_get( s->x, 0 ),
							gsl_vector_get( s->x, 1 ),
							gsl_vector_get( s->x, 2 ),
							gsl_blas_dnrm2( s->f ) );
				}
				else
				{
					if ( profile_ == "EMG" )
					{
						printf ( "before loop iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
								gsl_vector_get( s->x, 0 ),
								gsl_vector_get( s->x, 1 ),
								gsl_vector_get( s->x, 2 ),
								gsl_vector_get( s->x, 3 ),
								gsl_blas_dnrm2( s->f ) );
					}
					else
					{
						printf ( "before loop iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
								gsl_vector_get( s->x, 0 ),
								gsl_vector_get( s->x, 1 ),
								gsl_vector_get( s->x, 2 ),
								gsl_vector_get( s->x, 3 ),
								//gsl_vector_get(s->x,4),
								gsl_blas_dnrm2( s->f ) );

					}
				}
#endif

				// this is the loop for fitting
				do
				{
					iter++;
					status = gsl_multifit_fdfsolver_iterate ( s );

#ifdef DEBUG_FEATUREFINDER
					// This is useful for debugging
					if ( profile_ == "LmaGauss" )
					{
						printf ( "in loop iter: %4u x = % 15.8f % 15.8f % 15.8f |f(x)| = %g\n", iter,
								gsl_vector_get( s->x, 0 ),
								gsl_vector_get( s->x, 1 ),
								gsl_vector_get( s->x, 2 ),
								gsl_blas_dnrm2( s->f ) );
					}
					else
					{
						if ( profile_ == "EMG" )
						{
							printf ( "in loop iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
									gsl_vector_get( s->x, 0 ),
									gsl_vector_get( s->x, 1 ),
									gsl_vector_get( s->x, 2 ),
									gsl_vector_get( s->x, 3 ),
									gsl_blas_dnrm2( s->f ) );
						}
						else
						{
							printf ( "in loop iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
									gsl_vector_get( s->x, 0 ),
									gsl_vector_get( s->x, 1 ),
									gsl_vector_get( s->x, 2 ),
									gsl_vector_get( s->x, 3 ),
									//gsl_vector_get(s->x,4),
									gsl_blas_dnrm2( s->f ) );
						}
					}
#endif

					// fit is done
					if ( status ) break;
					status = gsl_multifit_test_delta( s->dx, s->x, eps_abs_, eps_rel_ );
				}
				while ( status == GSL_CONTINUE && iter < max_iteration_ );

				// This function uses Jacobian matrix J to compute the covariance matrix of the best-fit parameters, covar. The parameter epsrel (0.0) is used to remove linear-dependent columns when J is rank deficient.
				gsl_multifit_covar( s->J, 0.0, covar );

#ifdef DEBUG_FEATUREFINDER
				gsl_matrix_fprintf( stdout, covar, "%g" );
#endif

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

				gsl_status_ = gsl_strerror( status );

#ifdef DEBUG_FEATUREFINDER
				std::cout << profile_ << " status: " << gsl_status_ << std::endl;
#endif

				if ( profile_ == "LmaGauss" )
				{
#ifdef DEBUG_FEATUREFINDER
					printf( "deviation          = %.5f +/- %.5f\n", FIT( 0 ), ERR( 0 ) );
					printf( "expected_value	   = %.5f +/- %.5f\n", FIT( 1 ), ERR( 1 ) );
					printf( "scale_factor       = %.5f +/- %.5f\n", FIT( 2 ), ERR( 2 ) );
#endif
					standard_deviation_ = FIT( 0 );
					expected_value_ = FIT( 1 );
					scale_factor_ = FIT( 2 );
				}
				else
				{
					if ( profile_ == "EMG" )
					{
#ifdef DEBUG_FEATUREFINDER
						printf( "h = %.5f +/- %.5f\n", FIT( 0 ), ERR( 0 ) );
						printf( "w = %.5f +/- %.5f\n", FIT( 1 ), ERR( 1 ) );
						printf( "s = %.5f +/- %.5f\n", FIT( 2 ), ERR( 2 ) );
						printf( "z = %.5f +/- %.5f\n", FIT( 3 ), ERR( 3 ) );
#endif
						height_ = FIT( 0 );
						width_ = FIT( 1 );
						symmetry_ = FIT( 2 );
						retention_ = FIT( 3 );
					}
					else
					{
#ifdef DEBUG_FEATUREFINDER
						printf( "h = %.5f +/- %.5f\n", FIT( 0 ), ERR( 0 ) );
						printf( "w = %.5f +/- %.5f\n", FIT( 1 ), ERR( 1 ) );
						printf( "s = %.5f +/- %.5f\n", FIT( 2 ), ERR( 2 ) );
						printf( "z = %.5f +/- %.5f\n", FIT( 3 ), ERR( 3 ) );
						//	printf("r = %.5f +/- %.5f\n", FIT(4), ERR(4));
#endif
						height_ = FIT( 0 );
						width_ = FIT( 1 );
						symmetry_ = FIT( 2 );
						retention_ = FIT( 3 );
						r_ = r_; //FIT(4);
					}
				}

#ifdef DEBUG_FEATUREFINDER

				{
					// chi-squared value
					double chi = gsl_blas_dnrm2( s->f );
					printf( "chisq/dof = %g\n", pow( chi, 2.0 ) / ( n - p ) );
				}
#endif

				// function free all memory associated with the solver s
				gsl_multifit_fdfsolver_free( s );

#ifdef DEBUG_FEATUREFINDER
				for ( size_t current_point = 0; current_point < positionsDC_.size();current_point++ )
					std::cout << positionsDC_[ current_point ] << " " << signalDC_[ current_point ] << std::endl;

				std::cout << "" << std::endl;
				std::cout << "*** parameter for optimization ***" << std::endl;
				std::cout << "       height:  " << height_ << std::endl;
				std::cout << "        width:  " << width_ << std::endl;
				std::cout << "     symmetry:  " << symmetry_ << std::endl;
				std::cout << "    retention:  " << retention_ << std::endl;
				std::cout << "std.deviation:  " << standard_deviation_ << std::endl;
				std::cout << "max_iteration:  " << max_iteration_ << std::endl;
				std::cout << "      eps_abs:  " << eps_abs_ << std::endl;
				std::cout << "      eps_rel:  " << eps_rel_ << std::endl;
				std::cout << "      profile:  " << profile_ << std::endl;
				std::cout << "" << std::endl;
#endif
				positionsDC_.clear();
				signalDC_.clear();
			}

			/// get height for the EMG and logNormal model
			CoordinateType getHeight() const
			{
				return height_;
			}
			/// get width for the EMG and logNormal model
			CoordinateType getWidth() const
			{
				return width_;
			}
			/// get symmetry for the EMG and logNormal model
			CoordinateType getSymmetry() const
			{
				return symmetry_;
			}
			/// get retention time for the EMG and logNormal model
			CoordinateType getRT() const
			{
				return retention_;
			}
			/// get standard deviation for the Gauss
			CoordinateType getStandardDeviation() const
			{
				return standard_deviation_;
			}
			/// get expected value for the Gauss
			CoordinateType getExpectedValue() const
			{
				return expected_value_;
			}
			/// get scale factor for the Gauss
			CoordinateType getScaleFactor() const
			{
				return scale_factor_;
			}
			/// get GSL status
			std::string getGSLStatus() const
			{
				return gsl_status_;
			}

			/// Evaluation of the target function for nonlinear optimization.
			static int residualDC(const gsl_vector* x, void* params, gsl_vector* f)
			{
				size_t n = ( ( struct Internal::ExpFitPolyData* ) params ) ->n;
				String profile = ( ( struct Internal::ExpFitPolyData* ) params ) ->profile;

				/// normal distribution (s = standard deviation, m = expected value)
				if ( profile == "LmaGauss" )
				{
					double normal_s = gsl_vector_get( x, 0 );
					double normal_m = gsl_vector_get( x, 1 );
					double normal_scale = gsl_vector_get( x, 2 );

					double Yi = 0.0;

					for ( size_t i = 0; i < n; i++ )
					{
						double t = positionsDC_[ i ];

						Yi = ( 1 / ( sqrt( 2 * M_PI ) * normal_s ) ) * exp( -( ( t - normal_m ) * ( t - normal_m ) ) / ( 2 * normal_s * normal_s ) ) * normal_scale;

						gsl_vector_set( f, i, ( Yi - signalDC_[ i ] ) );
					}
				}
				else
				{
					/// Simplified EMG
					if ( profile == "EMG" )
					{
						double h = gsl_vector_get( x, 0 );
						double w = gsl_vector_get( x, 1 );
						double s = gsl_vector_get( x, 2 );
						double z = gsl_vector_get( x, 3 );

						double Yi = 0.0;

						// iterate over all points of the signal
						for ( size_t i = 0; i < n; i++ )
						{
							double t = positionsDC_[ i ];

							// Simplified EMG
							Yi = ( h * w / s ) * sqrt( 2 * M_PI ) * exp( ( pow( w, 2 ) / ( 2 * pow( s, 2 ) ) ) - ( ( t - z ) / s ) ) / ( 1 + exp( ( -2.4055 / sqrt( 2 ) ) * ( ( ( t - z ) / w ) - w / s ) ) );

							gsl_vector_set( f, i, ( Yi - signalDC_[ i ] ) );
						}
					}
					/// log normal
					else
					{
						double h = gsl_vector_get( x, 0 );
						double w = gsl_vector_get( x, 1 );
						double s = gsl_vector_get( x, 2 );
						double z = gsl_vector_get( x, 3 );
						double r = 2; //gsl_vector_get(x,4);

						double Yi = 0.0;

						for ( size_t i = 0; i < n; i++ )
						{
							double t = positionsDC_[ i ];

							Yi = h * exp( -log( r ) / ( log( s ) * log( s ) ) * pow( log( ( t - z ) * ( s * s - 1 ) / ( w * s ) + 1 ), 2 ) );

							gsl_vector_set( f, i, ( Yi - signalDC_[ i ] ) );
						}
					}
				}

				return GSL_SUCCESS;
			}

			/// Compute the Jacobian of the residual, where each row of the matrix corresponds to a point in the data.
			static int jacobianDC(const gsl_vector* x, void* params, gsl_matrix* J)
			{

				size_t n = ( ( struct Internal::ExpFitPolyData* ) params ) ->n;
				String profile = ( ( struct Internal::ExpFitPolyData* ) params ) ->profile;

				// normal distribution (s = standard deviation, m = expected value)
				if ( profile == "LmaGauss" )
				{
					double normal_s = gsl_vector_get( x, 0 );
					double normal_m = gsl_vector_get( x, 1 );
					double normal_scale = gsl_vector_get( x, 2 );

					double derivative_normal_s, derivative_normal_m, derivative_normal_scale = 0.0;

					for ( size_t i = 0; i < n; i++ )
					{
						double t = positionsDC_[ i ];

						// f'(normal_s)
						derivative_normal_s = -( ( 1 / sqrt( 2 * M_PI ) ) / ( normal_s * normal_s ) ) * exp( -( ( t - normal_m ) * ( t - normal_m ) ) / ( 2 * normal_s * normal_s ) ) * normal_scale + ( ( 1 / sqrt( 2 * M_PI ) ) / ( normal_s * normal_s * normal_s * normal_s ) ) * ( ( t - normal_m ) * ( t - normal_m ) ) * exp( -( ( t - normal_m ) * ( t - normal_m ) ) / ( 2 * normal_s * normal_s ) ) * normal_scale;

						// f'(normal_m)
						derivative_normal_m = ( ( 1 / sqrt( 2 * M_PI ) ) / ( normal_s * normal_s * normal_s ) ) * ( t - normal_m ) * exp( -( ( t - normal_m ) * ( t - normal_m ) ) / ( 2 * normal_s * normal_s ) ) * normal_scale;

						// f'(normal_scale)
						derivative_normal_scale = ( ( 1 / sqrt( 2 * M_PI ) ) / ( normal_s ) ) * exp( -( ( t - normal_m ) * ( t - normal_m ) ) / ( 2 * normal_s * normal_s ) );

						// set the jacobian matrix of the normal distribution
						gsl_matrix_set( J, i, 0, derivative_normal_s );
						gsl_matrix_set( J, i, 1, derivative_normal_m );
						gsl_matrix_set( J, i, 2, derivative_normal_scale );
					}
				}
				else
				{
					//Simplified EMG (sEMG)
					if ( profile == "EMG" )
					{
						double h = gsl_vector_get( x, 0 );
						double w = gsl_vector_get( x, 1 );
						double s = gsl_vector_get( x, 2 );
						double z = gsl_vector_get( x, 3 );

						const double emg_const = 2.4055;
						const double sqrt_2pi = sqrt( 2 * M_PI );
						const double sqrt_2 = sqrt( 2 );

						double exp1, exp2, exp3 = 0.0;
						double derivative_height, derivative_width, derivative_symmetry, derivative_retention = 0.0;

						// iterate over all points of the signal
						for ( size_t i = 0; i < n; i++ )
						{
							double t = positionsDC_[ i ];

							exp1 = exp( ( ( w * w ) / ( 2 * s * s ) ) - ( ( t - z ) / s ) );
							exp2 = ( 1 + exp( ( -emg_const / sqrt_2 ) * ( ( ( t - z ) / w ) - w / s ) ) );
							exp3 = exp( ( -emg_const / sqrt_2 ) * ( ( ( t - z ) / w ) - w / s ) );

							// f'(h) - sEMG
							derivative_height = w / s * sqrt_2pi * exp1 / exp2;

							// f'(h) - sEMG
							derivative_width = h / s * sqrt_2pi * exp1 / exp2 + ( h * w * w ) / ( s * s * s ) * sqrt_2pi * exp1 / exp2 + ( emg_const * h * w ) / s * sqrt_2pi * exp1 * ( -( t - z ) / ( w * w ) - 1 / s ) * exp3 / ( ( exp2 * exp2 ) * sqrt_2 );

							// f'(s) - sEMG
							derivative_symmetry = - h * w / ( s * s ) * sqrt_2pi * exp1 / exp2 + h * w / s * sqrt_2pi * ( -( w * w ) / ( s * s * s ) + ( t - z ) / ( s * s ) ) * exp1 / exp2 + ( emg_const * h * w * w ) / ( s * s * s ) * sqrt_2pi * exp1 * exp3 / ( ( exp2 * exp2 ) * sqrt_2 );

							// f'(z) - sEMG
							derivative_retention = h * w / ( s * s ) * sqrt_2pi * exp1 / exp2 - ( emg_const * h ) / s * sqrt_2pi * exp1 * exp3 / ( ( exp2 * exp2 ) * sqrt_2 );

							// set the jacobian matrix
							gsl_matrix_set( J, i, 0, derivative_height );
							gsl_matrix_set( J, i, 1, derivative_width );
							gsl_matrix_set( J, i, 2, derivative_symmetry );
							gsl_matrix_set( J, i, 3, derivative_retention );
						}
					}
					/// log normal function
					else
					{
						double h = gsl_vector_get( x, 0 );
						double w = gsl_vector_get( x, 1 );
						double s = gsl_vector_get( x, 2 );
						double z = gsl_vector_get( x, 3 );
						double r = 2; //gsl_vector_get(x,4);

						double derivative_height, derivative_width, derivative_symmetry, derivative_retention, derivative_r = 0.0;

						// iterate over all points of the signal
						for ( size_t i = 0; i < n; i++ )
						{
							double t = positionsDC_[ i ];

							double exp1 = exp( -log( r ) / ( log( s ) * log( s ) ) * pow( log( ( t - z ) * ( s * s - 1 ) / ( w * s ) + 1 ), 2 ) );
							double term1 = ( ( ( t - z ) * ( s * s - 1 ) ) / ( w * s ) ) + 1;
							double log_s = log( s );
							double log_term1 = log( term1 );
							double log_r = log( r );

							derivative_height = exp1;

							derivative_width = 2 * h * log_r / ( log_s * log_s ) * log_term1 * ( t - z ) * ( s * s - 1 ) / ( w * w ) / s / term1 * exp1;

							derivative_symmetry = h * ( 2 * log_r / ( log_s * log_s * log_s ) * ( log_term1 * log_term1 ) / s - 2 * log_r / ( log_s * log_s ) * log_term1 * ( 2 * ( t - z ) / w - ( t - z ) * ( s * s - 1 ) / ( w * s * s ) ) / term1 ) * exp1;

							derivative_retention = 2 * h * log_r / ( log_s * log_s ) * log_term1 * ( s * s - 1 ) / ( w * s ) / term1 * exp1;

							derivative_r = -h / r / ( log_s * log_s ) * ( log_term1 * log_term1 ) * exp1;

							// set the jacobian matrix
							gsl_matrix_set( J, i, 0, derivative_height );
							gsl_matrix_set( J, i, 1, derivative_width );
							gsl_matrix_set( J, i, 2, derivative_symmetry );
							gsl_matrix_set( J, i, 3, derivative_retention );
							//gsl_matrix_set(J, i, 4, derivative_r);
						}
					}
				}

				return GSL_SUCCESS;
			} // jacobianDC() 

			/// Driver function for the evaluation of function and jacobian.
			static int evaluateDC(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
			{
				residualDC( x, params, f );
				jacobianDC( x, params, J );

				return GSL_SUCCESS;
			}



		protected:

			virtual void updateMembers_()
			{
				max_iteration_ = this->param_.getValue( "rt:max_iteration" );
				eps_abs_ = this->param_.getValue( "rt:deltaAbsError" );
				eps_rel_ = this->param_.getValue( "rt:deltaRelError" );
				profile_ = ( std::string ) this->param_.getValue( "rt:profile" );

				interpolation_step_mz_ = this->param_.getValue( "mz:interpolation_step" );
				interpolation_step_rt_ = this->param_.getValue( "rt:interpolation_step" );

				iso_stdev_first_ = this->param_.getValue( "isotope_model:stdev:first" );
				iso_stdev_last_ = this->param_.getValue( "isotope_model:stdev:last" );
				iso_stdev_stepsize_ = this->param_.getValue( "isotope_model:stdev:step" );

				first_mz_model_ = ( Int ) this->param_.getValue( "mz:model_type:first" );
				last_mz_model_ = ( Int ) this->param_.getValue( "mz:model_type:last" );
			}


			/// fit offset by maximizing of quality
			double fitOffset_(InterpolationModel* model, const IndexSet& set, double stdev1, double stdev2, Coordinate offset_step)
			{
				const Coordinate offset_min = model->getInterpolation().supportMin() - stdev1;
				const Coordinate offset_max = model->getInterpolation().supportMin() + stdev2;

				Coordinate offset;
				QualityType correlation;

				//test model with default offset
				std::vector<Real> data;
				data.reserve(set.size());
				std::vector<Real> model_data;
				model_data.reserve(set.size());
				for (IndexSet::iterator it=set.begin();it!=set.end();++it)
				{
					data.push_back(this->getPeakIntensity(*it));
					model_data.push_back(model2D_.getIntensity(DPosition<2>(this->getPeakRt(*it),this->getPeakMz(*it))));
				}

				Coordinate max_offset = model->getInterpolation().getOffset();
				QualityType max_correlation = Math::BasicStatistics<Real>::pearsonCorrelationCoefficient(data.begin(), data.end(), model_data.begin(), model_data.end());

				//test different offsets
				for ( offset = offset_min; offset <= offset_max; offset += offset_step )
				{
					model->setOffset( offset );
					model_data.clear();
					for (IndexSet::iterator it=set.begin();it!=set.end();++it)
					{
						model_data.push_back(model2D_.getIntensity(DPosition<2>(this->getPeakRt(*it),this->getPeakMz(*it))));
					}
					// IWASHEREMSG(std::distance(data.begin(),data.end()));
					// IWASHEREMSG(std::distance(model_data.begin(),model_data.end()));
					correlation = Math::BasicStatistics<Real>::pearsonCorrelationCoefficient(data.begin(), data.end(), model_data.begin(), model_data.end());
					if ( correlation > max_correlation )
					{
						max_correlation = correlation;
						max_offset = offset;
					}
				}
				model->setOffset( max_offset );
				return max_correlation;
			}


			double fit_(const IndexSet& set, MzFitting mz_fit, RtFitting rt_fit, Coordinate isotope_stdev=0.1)
			{
				// IWASHEREMSG(isotope_stdev);

				// Build Models
				InterpolationModel * mz_model;
				if ( mz_fit == MZGAUSS )
				{
					mz_model = new GaussModel();
					mz_model->setInterpolationStep( interpolation_step_mz_ );

					Param tmp;
					tmp.setValue( "bounding_box:min", min_[ MZ ] );
					tmp.setValue( "bounding_box:max", max_[ MZ ] );
					tmp.setValue( "statistics:variance", mz_stat_.variance() );
					tmp.setValue( "statistics:mean", mz_stat_.mean() );

					static_cast<GaussModel*>( mz_model ) ->setParameters( tmp );
				}
				else
				{
					// new model
					mz_model = new IsotopeModel();
					Param iso_param = this->param_.copy( "isotope_model:", true );
					iso_param.remove( "stdev" );
					mz_model->setParameters( iso_param );
					mz_model->setInterpolationStep( interpolation_step_mz_ );

					Param tmp;
					tmp.setValue( "charge", static_cast<Int>( mz_fit ) );
					tmp.setValue( "isotope:stdev", isotope_stdev );
					tmp.setValue( "statistics:mean", mz_stat_.mean() );

					static_cast<IsotopeModel*>( mz_model ) ->setParameters( tmp );
				}

				InterpolationModel* rt_model;
				if ( rt_fit == RTGAUSS )
				{
					rt_model = new GaussModel();
					rt_model->setInterpolationStep( interpolation_step_rt_ );

					Param tmp;
					tmp.setValue( "bounding_box:min", min_[ RT ] );
					tmp.setValue( "bounding_box:max", max_[ RT ] );
					tmp.setValue( "statistics:variance", rt_stat_.variance() );
					tmp.setValue( "statistics:mean", rt_stat_.mean() );

					static_cast<GaussModel*>( rt_model ) ->setParameters( tmp );
				}
				else if ( rt_fit == LMAGAUSS )
				{
					rt_model = new LmaGaussModel();
					rt_model->setInterpolationStep( interpolation_step_rt_ );

					Param tmp;
					tmp.setValue( "bounding_box:min", min_[ RT ] );
					tmp.setValue( "bounding_box:max", max_[ RT ] );
					tmp.setValue( "statistics:variance", rt_stat_.variance() );
					tmp.setValue( "statistics:mean", rt_stat_.mean() );
					tmp.setValue( "lma:scale_factor", scale_factor_ );
					tmp.setValue( "lma:standard_deviation", standard_deviation_ );
					tmp.setValue( "lma:expected_value", expected_value_ );

					static_cast<LmaGaussModel*>( rt_model ) ->setParameters( tmp );
				}
				else if ( rt_fit == EMGAUSS )
				{
					rt_model = new EmgModel();
					rt_model->setInterpolationStep( interpolation_step_rt_ );

					Param tmp;
					tmp.setValue( "bounding_box:min", min_[ RT ] );
					tmp.setValue( "bounding_box:max", max_[ RT ] );
					tmp.setValue( "statistics:variance", rt_stat_.variance() );
					tmp.setValue( "statistics:mean", rt_stat_.mean() );
					tmp.setValue( "emg:height", height_ );
					tmp.setValue( "emg:width", width_ );
					tmp.setValue( "emg:symmetry", symmetry_ );
					tmp.setValue( "emg:retention", retention_ );

					static_cast<LmaGaussModel*>( rt_model ) ->setParameters( tmp );
				}
				else if ( rt_fit == LOGNORMAL )
				{
					rt_model = new LogNormalModel();
					rt_model->setInterpolationStep( interpolation_step_rt_ );

					Param tmp;
					tmp.setValue( "bounding_box:min", min_[ RT ] );
					tmp.setValue( "bounding_box:max", max_[ RT ] );
					tmp.setValue( "statistics:variance", rt_stat_.variance() );
					tmp.setValue( "statistics:mean", rt_stat_.mean() );
					tmp.setValue( "emg:height", height_ );
					tmp.setValue( "emg:width", width_ );
					tmp.setValue( "emg:symmetry", symmetry_ );
					tmp.setValue( "emg:retention", retention_ );
					tmp.setValue( "lognormal:r", r_ );

					static_cast<LmaGaussModel*>( rt_model ) ->setParameters( tmp );
				}
				else
				{
					rt_model = new BiGaussModel();
					rt_model->setInterpolationStep( interpolation_step_rt_ );

					Param tmp;
					tmp.setValue( "bounding_box:min", min_[ RT ] );
					tmp.setValue( "bounding_box:max", max_[ RT ] );
					tmp.setValue( "statistics:mean", rt_stat_.mean() );
					tmp.setValue( "statistics:variance1", rt_stat_.variance1() );
					tmp.setValue( "statistics:variance2", rt_stat_.variance2() );

					static_cast<BiGaussModel*>( rt_model ) ->setParameters( tmp );
				}

				model2D_.setModel( MZ, mz_model ).setModel( RT, rt_model );

				QualityType res;
				StopWatch w;
				w.start();
				res = fitOffset_( mz_model, set, stdev_mz_, stdev_mz_, interpolation_step_mz_ );
				w.stop();
				std::cout << "Time spent for mz offset: " << w.getClockTime() << std::endl;

				if ( profile_ != "LmaGauss" && profile_ != "EMG" && profile_ != "LogNormal" )
				{
					res = fitOffset_( rt_model, set, stdev_rt1_, stdev_rt2_, interpolation_step_rt_ );
				}
				else
				{
					//???? debugging
					std::cerr << "Unrecognized profile: '" << profile_ << "'" << std::endl;	
				}
				return res;
			}

			ProductModel<2> model2D_;
			Math::BasicStatistics<> mz_stat_;
			Math::AsymmetricStatistics<> rt_stat_;
			double stdev_mz_;
			double stdev_rt1_;
			double stdev_rt2_;
			DPosition<2> min_;
			DPosition<2> max_;

			/// counts features (used for debug output only)
			UInt counter_;

			/// interpolation step size (in m/z)
			Coordinate interpolation_step_mz_;
			/// interpolation step size (in retention time)
			Coordinate interpolation_step_rt_;

			/// first stdev
			float iso_stdev_first_;
			/// last stdev
			float iso_stdev_last_;
			/// step size
			float iso_stdev_stepsize_;

			/// first mz model (0=Gaussian, 1....n = charge )
			Int first_mz_model_;			
			/// last mz model
			Int last_mz_model_;

			/// Maximum number of iterations
			unsigned int max_iteration_;

			/// parameter of log normal function:
			/// r is the ratio between h and the height at which w and s are computed
			double r_;

			/// parameter of emg and log normal function:height
			double height_;
			/// parameter of emg and log normal function: width
			double width_;
			/// parameter of emg and log normal function: symmetry
			double symmetry_;
			/// parameter of emg and log normal function: retention time
			double retention_;
			/// parameter indicates symmetric peaks
			bool symmetric_;
			/// gsl status
			std::string gsl_status_;
			/// function for fitting
			std::string profile_;

			/** Test for the convergence of the sequence by comparing the last iteration step dx with the absolute error epsabs and relative error epsrel to the current position x */
			/// absolute error
			double eps_abs_;
			/// relative error
			double eps_rel_;

			/// parameter of gauss function: standard deviation
			double standard_deviation_;
			/// parameter of gauss function: scale factor
			double scale_factor_;
			/// parameter of gauss function: expected value
			double expected_value_;		

			//positions and signal values
			static std::vector<double> positionsDC_;
			static std::vector<double> signalDC_;


		private:
			/// Not implemented
			SimpleModelFitter();
			/// Not implemented
			SimpleModelFitter& operator=(const SimpleModelFitter&);
			/// Not implemented
			SimpleModelFitter(const SimpleModelFitter&);

	};

	template<typename P,typename F>	std::vector<double> SimpleModelFitter<P,F>::positionsDC_;
	template<typename P,typename F>	std::vector<double> SimpleModelFitter<P,F>::signalDC_;
	
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEMODELFITTER_H
