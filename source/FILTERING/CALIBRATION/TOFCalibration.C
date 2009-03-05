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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/CALIBRATION/TOFCalibration.h>

namespace OpenMS
{
	
  TOFCalibration::TOFCalibration()
    :DefaultParamHandler("TOFCalibration"),ProgressLogger()
  {
    subsections_.push_back("PeakPicker");
		check_defaults_=false; // class has no own parameters
  }

  TOFCalibration::~TOFCalibration(){}


  void TOFCalibration::calculateCalibCoeffs_(MSExperiment<>& calib_spectra)
  {
		// flight times are needed later
    calib_peaks_ft_ = calib_spectra;


    // convert flight times of peaks into m/z values
    applyTOFConversion_(calib_spectra);
		std::vector<std::vector<unsigned int> > monoiso_peaks;
    getMonoisotopicPeaks_(calib_spectra,monoiso_peaks);
		
    startProgress(0,calib_spectra.size(),"quadratic fitting of calibrant spectra");
    // do the quadratic fitting for each calibration spectra separately
    for(unsigned int spec=0;spec<calib_spectra.size();++spec)
      {
				std::vector<unsigned int> monoiso_peaks_scan;
				std::vector<double> exp_masses;
				// match the m/z-values to the expected masses
				matchMasses_(calib_spectra,monoiso_peaks,monoiso_peaks_scan,exp_masses,spec);
				
				// the actual quadratic fitting part
				gsl_matrix *X, *cov;
				gsl_vector *y, *c;
				Size n = exp_masses.size();
				if(n < 3)
					{
						continue;
					}
				
				double chisq;
				
				// matrix containing the observations
				X = gsl_matrix_alloc (n, 3);
				// vector containing the expected masses
				y = gsl_vector_alloc (n);
				
				// vector containing the coefficients of the quadratic function after the fitting
				c = gsl_vector_alloc (3);
				// matrix containing the covariances
				cov = gsl_matrix_alloc (3, 3);
				
				for (Size i = 0; i < n; i++)
					{
						// get the flight time
						double xi = ( (calib_peaks_ft_.begin() + spec)->begin() + monoiso_peaks_scan[i])->getMZ();
						// y_i = a + b*x_i + c*x_i^2  <---- the quadratic equation, a, b, and c shall be determined
						// x_i^0 = 1 --> enter 1 at the first position of each row 
						gsl_matrix_set (X, i, 0, 1.0);

						// x_i^1 at the second position
						gsl_matrix_set (X, i, 1, xi);
						// x_i^2 at the third position
						gsl_matrix_set (X, i, 2, xi*xi);
						
						// set expected mass
						gsl_vector_set (y, i, exp_masses[i]);
						
					}
				
				gsl_multifit_linear_workspace * work
					= gsl_multifit_linear_alloc (n, 3);
				
				gsl_multifit_linear (X, y, c, cov,
														 &chisq, work);


#ifdef DEBUG_CALIBRATION
				printf ("# best fit: Y = %g + %g X + %g X^2\n",
								gsl_vector_get(c,(0)),
								gsl_vector_get(c,(1)),
								gsl_vector_get(c,(2)));
#endif
				// store the coefficients
				coeff_quad_fit_.push_back(gsl_vector_get(c,(0)));
				coeff_quad_fit_.push_back(gsl_vector_get(c,(1)));
				coeff_quad_fit_.push_back(gsl_vector_get(c,(2)));
				
				gsl_multifit_linear_free (work);

				// determine the errors in ppm

				
				std::vector<double> errors;
				for(Size p = 0; p < n; ++p)
					{
#ifdef DEBUG_CALIBRATION
						std::cout << exp_masses[p]
											<< "\t"<<mQ_(calib_peaks_ft_[spec][monoiso_peaks_scan[p]].getMZ(),spec) - exp_masses[p]<< std::endl;
#endif
						errors_[exp_masses[p]].push_back((mQ_(calib_peaks_ft_[spec][monoiso_peaks_scan[p]].getMZ(),spec) - exp_masses[p]));
					}
				setProgress(spec);
      }
    endProgress();

		if(coeff_quad_fit_.size() == 0)
			{
				String mess = String("Data can't be calibrated, not enough reference masses found: ") + coeff_quad_fit_.size()/3;
				Exception::UnableToCalibrate(__FILE__, __LINE__,__PRETTY_FUNCTION__,"UnableToCalibrate", mess.c_str());
			}
    averageErrors_();
    averageCoefficients_();
		
		
    double* calib_masses = new double[error_medians_.size()];
    double* error_medians = new double[error_medians_.size()];
    for(unsigned int i=0; i < error_medians_.size(); ++i )
      {
				calib_masses[i] = calib_masses_[i];
				error_medians[i] = error_medians_[i];				
      }

    acc_ = gsl_interp_accel_alloc();
    spline_ = gsl_spline_alloc(gsl_interp_cspline,error_medians_.size());
    gsl_spline_init(spline_,calib_masses,error_medians,error_medians_.size());
	
#ifdef DEBUG_CALIBRATION
    std::cout<< "fehler nach spline fitting" << std::endl;
		
    for(unsigned int spec=0;spec <  calib_peaks_ft_.size(); ++spec)
      {
				
				std::vector<double> exp_masses;
				std::vector<unsigned int> monoiso;
				matchMasses_(calib_spectra,monoiso_peaks,monoiso,exp_masses,spec);
				for(unsigned int p=0; p < monoiso.size();++p)
					{
						double xi = mQ_(calib_peaks_ft_[spec][monoiso[p]].getMZ(),spec);
						if(xi > calib_masses[error_medians_.size()-1] ) continue;
						if(xi < calib_masses[0]) continue;
						std::cout << exp_masses[p] << "\t" 
											<< (xi - exp_masses[p] - gsl_spline_eval(spline_,xi,acc_))/exp_masses[p] * 1e6
											<< std::endl;
						
					}

      }


    double xi,yi;
    std::cout << "interpolation \n\n";
    for(xi = calib_masses[0]; xi < calib_masses[error_medians_.size()-1]; xi += 0.01)
      {
				yi = gsl_spline_eval(spline_,xi,acc_);
				std::cout << xi << "\t" << yi << std::endl;
      }
    std::cout << "--------------\nend interpolation \n\n";
#endif

    delete[] calib_masses;
    delete[] error_medians;
  }


	

	
  void TOFCalibration::averageCoefficients_()
  {
    a_ = 0;
    b_ = 0;
    c_ = 0;
    for(unsigned int i=0;i < coeff_quad_fit_.size(); i+=3)
      {
				a_ += coeff_quad_fit_[i];
				b_ += coeff_quad_fit_[i+1];
				c_ += coeff_quad_fit_[i+2];
      }
    a_/=(coeff_quad_fit_.size()/3);
    b_/=(coeff_quad_fit_.size()/3);
    c_/=(coeff_quad_fit_.size()/3);
  }
	
  void TOFCalibration::averageErrors_()
  {
    for(unsigned int p=0;p<exp_masses_.size();++p)
      {
				// mean
				if(errors_[exp_masses_[p]].size()>0)
					{
						double sum=0;
						for(unsigned int i = 0; i < errors_[exp_masses_[p]].size() ; ++i)
							{
								sum += errors_[exp_masses_[p]][i];
							
							}
						error_medians_.push_back(sum/errors_[exp_masses_[p]].size());
						calib_masses_.push_back(exp_masses_[p]);
					}
      }
  }

	

  void TOFCalibration::matchMasses_(MSExperiment<>& calib_peaks,
																				 std::vector<std::vector<unsigned int> >& monoiso_peaks,
																				 std::vector<unsigned int>& obs_masses,
																				 std::vector<double>& exp_masses,unsigned int idx)
  {
    for(unsigned int i=0;i<monoiso_peaks[idx].size();++i)
      {
				for(unsigned int j=0; j < exp_masses_.size(); ++j)
					{
						if(fabs(( (calib_peaks.begin() + idx)->begin() + (monoiso_peaks[idx])[i])->getMZ() - exp_masses_[j] ) < 1)
							{
								obs_masses.push_back((monoiso_peaks[idx])[i]);
								exp_masses.push_back(exp_masses_[j]);
								break;
							}
					}
      }
#ifdef DEBUG_CALIBRATION

    std::cout << "\n\n---------\nmatching monoisotopic peaks\n";
		
    for(unsigned int i=0;i<obs_masses.size();++i)
      {
				std::cout << ( (calib_peaks_ft_.begin() + idx)->begin() + obs_masses[i])->getMZ()
									<< "\t"<<exp_masses[i]
									<< std::endl;
				
      }

#endif
  }


  void TOFCalibration::getMonoisotopicPeaks_(MSExperiment<>& calib_peaks, std::vector<std::vector<unsigned int> >& monoiso_peaks)
  {
		
    MSExperiment<>::iterator spec_iter = calib_peaks.begin();
    MSExperiment<>::SpectrumType::iterator peak_iter, help_iter;
#ifdef DEBUG_CALIBRATION
    spec_iter = calib_peaks.begin();
    std::cout << "\n\nbefore---------\n\n";
    // iterate through all spectra
    for(;spec_iter != calib_peaks.end();++spec_iter)
      {
				peak_iter = spec_iter->begin();
				// go through current scan
				for(;peak_iter != spec_iter->end();++peak_iter)
					{
						std::cout << peak_iter->getMZ() << std::endl;
					}
      }

#endif
    spec_iter = calib_peaks.begin();
    // iterate through all spectra
    for(;spec_iter != calib_peaks.end();++spec_iter)
      {
				peak_iter = spec_iter->begin();
				help_iter = peak_iter;
				std::vector<unsigned int> vec;
				// go through current scan
				while(peak_iter < spec_iter->end())
					{
						while(peak_iter+1 < spec_iter->end() && ( (peak_iter+1)->getMZ() - peak_iter->getMZ() < 1.2) )
							{
								++peak_iter;
							}
						
						vec.push_back(distance(spec_iter->begin(),help_iter));
					
						help_iter = peak_iter+1;
						++peak_iter;
						
					}
				monoiso_peaks.push_back(vec);

      }
		
#ifdef DEBUG_CALIBRATION

		
    std::cout << "\n\nafter---------\n\n";

    for(unsigned int i=0;i<monoiso_peaks.size();++i)
      {
				for(unsigned int j=0;j<monoiso_peaks[i].size();++j)
					{
						std::cout <<i<<"\t" <<( (calib_peaks.begin() + i)->begin() + (monoiso_peaks[i])[j])->getMZ() << std::endl;
					}
				std::cout << "--------------\n";
						
      }
    std::cout << "--------------\n\n\n";
#endif
  }

  void TOFCalibration::applyTOFConversion_(MSExperiment<>& calib_spectra)
  {
    MSExperiment<>::iterator spec_iter = calib_spectra.begin();
    MSExperiment<>::SpectrumType::iterator peak_iter;
    unsigned int idx =0;

		//two point conversion
		if(ml3s_.size()==0)
			{
				for(;spec_iter != calib_spectra.end();++spec_iter)
					{
						peak_iter = spec_iter->begin();
						double ml1,ml2;
						if(ml1s_.size()==1)
							{
								ml1 = ml1s_[0];
								ml2 = ml2s_[0];
							}
						else
							{
								ml1 = ml1s_[idx];
								ml2 = ml2s_[idx];
							}
						
						// go through current scan
						for(;peak_iter != spec_iter->end();++peak_iter)
							{
								double time = peak_iter->getMZ();
								peak_iter->setPos( ml1/1E12 *(time *1000 - ml2) );
							}
						++idx;
					}
			}
		else
			{
				// three point conversion
				for(;spec_iter != calib_spectra.end();++spec_iter)
					{
						peak_iter = spec_iter->begin();
						double ml1,ml2,ml3;
						if(ml1s_.size()==1)
							{
								ml1 = ml1s_[0];
								ml2 = ml2s_[0];
								ml3 = ml3s_[0];
							}
						else
							{
								ml1 = ml1s_[idx];
								ml2 = ml2s_[idx];
								ml3 = ml3s_[idx];
							}
						
						// go through current scan
						for(;peak_iter != spec_iter->end();++peak_iter)
							{
								double time = peak_iter->getMZ();
								peak_iter->setPos( (-ml2 - (0.1E7 * (-5E5 + sqrt(0.25E12 - ml1*ml2*ml3 + ml1*ml3*time) ))/(ml1*ml3) +time) / ml3);
							}
						++idx;
					}
			}
		
  }
	
} //namespace openms
