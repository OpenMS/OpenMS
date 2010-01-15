// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Clemens Groepl $
// --------------------------------------------------------------------------

#include <numeric>
#include <math.h>

#include <OpenMS/SIMULATION/ElutionModel.h>

namespace OpenMS
{
    ElutionModel::ElutionModel()
		: InterpolationModel()
    {
      setName(getProductName());

			// Since the interpolation table is (re-)initialized after setting
			// parameters, we set an empty bounding_box to avoid silly computations.
			defaults_.setValue("bounding_box:min",0.0f,"Lower end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
			defaults_.setValue("bounding_box:max",0.0f,"Upper end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
			defaults_.setValue("interpolation_step",0.1,"Sampling rate for the interpolation of the model function.", StringList::create("advanced"));

      defaults_.setValue("statistics:mean",0.0f,"Centroid position of the model.", StringList::create("advanced"));
      defaults_.setValue("statistics:variance",1.0f,"The variance of the model.", StringList::create("advanced"));
      defaults_.setValue("emg:height",0.0f,"Height of the exponentially modified Gaussian.", StringList::create("advanced"));
      defaults_.setValue("emg:width",1.0f,"Width of the exponentially modified Gaussian.", StringList::create("advanced"));
      defaults_.setValue("emg:symmetry",0.0f,"Symmetry of the exponentially modified Gaussian.", StringList::create("advanced"));
      defaults_.setValue("emg:retention",0.0f,"Retention time of the exponentially modified Gaussian.", StringList::create("advanced"));
      
      defaultsToParam_();
    }

    ElutionModel::ElutionModel(const ElutionModel& source)
    : InterpolationModel(source)
    {
      setParameters( source.getParameters() );
      updateMembers_();
    }
  
    ElutionModel::~ElutionModel()
    {	
    }
  
    ElutionModel& ElutionModel::operator = (const ElutionModel& source)
    {
      if (&source == this) return *this;
      
      InterpolationModel::operator = (source);
      setParameters( source.getParameters() );
      updateMembers_();
      
      return *this;
    }

    void ElutionModel::setSamples()
    {
			ContainerType& data = interpolation_.getData();
      data.clear();
      if (max_==min_) return;
      data.reserve( UInt ( (max_-min_) / interpolation_step_ + 1 ) );
      CoordinateType pos = min_;
			
			// write parameter values to std::cout


 #ifdef DEBUG_SIM
			static int invocation_counter = 0;
			++invocation_counter;

			std::cout << "emg_invocation: "
								<< retention_ <<' '
								<< width_ <<' '
								<< height_ <<' '
								<< symmetry_ <<' '
								<< min_ <<' '
								<< max_ <<' '
								<< std::endl;
#endif

			// EMG is like a Gaussian if skewness is small.  This case distinctions
			// avoids numerical instabilities until I find another formula that
			// works for all parameter settings.  Any suggestions are welcome!
			// groepl@inf.fu-berlin.de
			// TODO: This was int abs(int). Not intended, I suppose ?
			if ( fabs(symmetry_/width_) > 0.1 )
			{
				// Clemens' code:
				// This version of the exponentially modified Gaussian (EMG) function is from
				// http://www.systat.com/products/tablecurve2d/help/?sec=1132
				// In Gnuplot notation it reads:
				//
				//  emg(a,b,c,d,x)=(a*c*sqrt_pi_2/d)*exp((b-x)/d+(c*c)/(2*d*d)) * ( d/abs(d) - erf( ((b-x)/c+c/d)/sqrt_2 ) )
				//
				// IMPORTANT NOTE: We changed the semantics of the symmetry parameter
				// such that positive symmetry indicates tailed peaks and negative
				// symmetry indicates fronted peaks.  Hence the sign of symmetry is
				// reversed in the formulas below:

				// symmetry_ = -symmetry_;

				//
				// Amplitude, deconvolved Gaussian = a
				// Center, deconvolved Gaussian = b
				// Area, devonvolved Gaussian of Function = sqrt(2*pi)*a*c
				// FWHM, deconvolved Gaussian = 2*sqrt(2*ln(2))*c
				// Time Constant, Exponential = d
				// Constraints: c > 0 and d != 0
				// 
				// sqrt(PI/2) = 1.253314137315500251
				// 1/sqrt(2*PI) = 0.3989422804014326779
				// 1/sqrt(2) = 0.7071067811865475244
				// PI/2 = 1.570796326794896619
			
				// emg(a,b,c,d,x) = ...

				// a*c*sqrt_pi_2/d
				const DoubleReal factor1 = 1.253314137315500251 * height_ * width_ / symmetry_;
			
				// c/(sqrt(2)*d) 
				const DoubleReal summand = 0.7071067811865475244 * width_ / symmetry_;
			
				// (c*c)/(2*d*d)
				const DoubleReal summand_squared = 0.5 * (width_*width_) / (symmetry_*symmetry_);

				// d/abs(d)
				const DoubleReal sign_of_symmetry_ = ( symmetry_> 0 ? 1 : -1 );

				// 1/d
				const DoubleReal factor2 = 1.0 / symmetry_;

				// 1/(sqrt(2)*c)
				const DoubleReal factor3 = 0.7071067811865475244 / width_;
			
				pos = min_;
				for ( UInt i = 0; pos < max_; ++i)
				{
					pos = min_ + i * interpolation_step_;
					DoubleReal tmp = pos - retention_;

					// Clemens' code:
					// data.push_back (EMG as described above)
					const DoubleReal emg_value =
						(
						 factor1
						 *
						 exp(
								 tmp * factor2
								 + 
								 summand_squared
								)
						 *
						 (
							sign_of_symmetry_
							-
              boost::math::tr1::erf(
							//erf(
									tmp * factor3
									+
									summand
								 )
						 )
						);
					data.push_back(emg_value);
				}
			}
			else //  otherwise,  width_/symmetry  is very small
			{
				for ( UInt i = 0; pos< max_; ++i)
				{
					pos = min_ + i * interpolation_step_;
					DoubleReal tmp = pos - retention_;
					// actually use a Gaussian
					const DoubleReal gauss_value = height_ * exp( - 0.5 * tmp * tmp / width_ / width_ );
					data.push_back(gauss_value);
				}
			} // end of case distinction

//  #ifdef DEBUG_SIM
// 			{
// 				const DoubleReal factor1 = 1.253314137315500251 * height_ * width_ / symmetry_;
// 				const DoubleReal summand = 0.7071067811865475244 * width_ / symmetry_;
// 				const DoubleReal summand_squared = 0.5 * (width_*width_) / (symmetry_*symmetry_);
// 				const DoubleReal sign_of_symmetry_ = ( symmetry_> 0 ? 1 : -1 );
// 				const DoubleReal factor2 = 1.0 / symmetry_;
// 				const DoubleReal factor3 = 0.7071067811865475244 / width_;
// 				pos = min_;
// 				for ( UInt i = 0; pos< max_; ++i)
// 				{
// 					pos = min_ + i * interpolation_step_;
// 					DoubleReal tmp = pos - retention_;
// 					const DoubleReal emg_value =
// 						(factor1 * exp(tmp * factor2 + summand_squared) * (sign_of_symmetry_ - erf(tmp * factor3 + summand)));
// 					const DoubleReal gauss_value =
// 						height_ * exp( - 0.5 * tmp * tmp / width_ / width_ );
// 					std::cout << "emg_values: "
// 										<< invocation_counter <<' '
// 										<< tmp <<' '
// 										<< emg_value <<' '
// 										<< gauss_value
// 										<< std::endl;
// 				}
// 			}
// #endif

      interpolation_.setScale  ( interpolation_step_ );
      interpolation_.setOffset ( min_ );
			return;
		}

    void ElutionModel::setOffset(CoordinateType offset)
    {
      DoubleReal diff = offset - getInterpolation().getOffset();
      min_ += diff;
      max_ += diff;
      statistics_.setMean(statistics_.mean() + diff);
  
      InterpolationModel::setOffset(offset);
  
      param_.setValue("bounding_box:min", min_);
      param_.setValue("bounding_box:max", max_);
      param_.setValue("statistics:mean", statistics_.mean());
    }

    ElutionModel::CoordinateType ElutionModel::getCenter() const
    {
      return statistics_.mean();
    }
    
    void ElutionModel::updateMembers_()
    {
      InterpolationModel::updateMembers_();
  
      min_ = param_.getValue("bounding_box:min");
      max_ = param_.getValue("bounding_box:max");
      statistics_.setMean( param_.getValue("statistics:mean") );
      statistics_.setVariance(param_.getValue("statistics:variance"));
      height_ = param_.getValue("emg:height");
      width_ = param_.getValue("emg:width");
      symmetry_ = param_.getValue("emg:symmetry");
      retention_ = param_.getValue("emg:retention");
      
      setSamples();
    }
}


//////////////////////////////////////////////////////////////////////////////

// Hey, this place is for my private notes!  Do NOT read the stuff below!
// Clemens Groepl

#if 0 

//------------------------------------------------------------

// Marcel's code:
// DoubleReal sqrt_2pi = sqrt(2*M_PI);
// DoubleReal term_sq2 = (-2.4055/sqrt(2));
// DoubleReal part1    = (height_*width_/symmetry_);
// DoubleReal part2    = pow(width_,2)/(2*pow(symmetry_,2));
// DoubleReal part3    = width_/symmetry_;
//
// data.push_back (Simplified EMG)
// data.push_back((part1*sqrt_2pi*exp(part2-(tmp/symmetry_))/(1+exp(term_sq2*((tmp/width_)-part3)))));

//------------------------------------------------------------

// // The following calculation is taken from:
// //
// // Roland Delley: Series for the Exponentially Modified Gaussian Peak Shape.  Anal. Chem. 1985, 57, 388.
// //
// // But this does not work yet (maybe I got things wrong, maybe the paper)
// //
// // const DoubleReal St = width_ / symmetry_;
// // DoubleReal tmp = pos - retention_;
// const DoubleReal T = tmp / width_;
// const DoubleReal T_squared = T * T;
// // const DoubleReal St = width_ / symmetry_;
// const DoubleReal x = St - T;
// const DoubleReal x_squared = x * x;
// const DoubleReal fx = 
// 	St *
// 	(
// 	 1.253314137315500251
// 	 *
// 	 exp ( x_squared / 2.0 )
// 	 -
// 	 x *
// 	 ( 1 + ( x_squared / 3.0 ) *
// 		 ( 1 + ( x_squared / 5.0 ) *
// 			 ( 1 + ( x_squared / 7.0 ) *
// 				 ( 1 + ( x_squared / 9.0 ) *
// 					 ( 1 + ( x_squared / 11.0 )
// 					 )
// 				 )
// 			 )
// 		 )
// 	 )
// 	);
// // const DoubleReal YT = height_ * exp(-T_squared/2.0); 
// const DoubleReal YT = height_ * fx * exp(-T_squared/2.0); 
				
// std::cout << "emg: " << tmp <<' '<< emg_value <<' '<< YT <<' '<< x_squared <<' '<< T << std::endl;

//------------------------------------------------------------

#endif
