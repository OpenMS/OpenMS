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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <gsl/gsl_statistics.h>

using namespace std;

namespace OpenMS
{
	namespace Math
	{
		PosteriorErrorProbabilityModel::PosteriorErrorProbabilityModel()
		: DefaultParamHandler("PosteriorErrorProbabilityModel"), negative_prior_(0.5),max_gumbel_(0),max_gauss_(0),smallest_score_(0)
  	{
			defaults_.setValue("number_of_bins", 100, "Number of bins used for visualisation. Only needed if each iteration step of the EM-Algorithm will be visualized", StringList::create("advanced"));
			defaults_.setValue("output_plots","false","If true every step of the EM-algorithm will be written to a file as a gnuplot formula",StringList::create("advanced"));
			defaults_.setValidStrings("output_plots",StringList::create("true,false"));
			defaults_.setValue("output_name","", "if output_plots is on, the output files will be saved in the following manner: <output_name>scores.txt for the scores and <output_name>step_* for each step of the EM-algorithm e.g. output_name = /usr/home/OMSSA123_ then /usr/home/OMSSA123_scores.txt, /usr/home/OMSSA123_step_0, /usr/home/OMSSA123_step_1 ... will be written. If no directory is specified, e.g. instead of '/usr/home/OMSSA123_' just OMSSA123_, the files will be written into the working directory.",StringList::create("advanced,output file"));
			defaultsToParam_();			
		}
		
		PosteriorErrorProbabilityModel::~PosteriorErrorProbabilityModel()
		{
		}
				
		void PosteriorErrorProbabilityModel::fit( vector<double> & x_scores, vector<double>& probabilities)
		{	
			if(x_scores.empty())
			{
				return;
			}
			
			probabilities.clear();
			probabilities.resize(x_scores.size());
			//-------------------------------------------------------------
			// Initializing Parameters
			//-------------------------------------------------------------
			double minimum = x_scores[0];
			for(std::vector< double>::iterator it = x_scores.begin(); it < x_scores.end(); ++it)
			{
				if(minimum > *it)
				{
					minimum = *it;
				}
			}
			minimum = fabs(minimum);
			vector< double>::iterator probs = probabilities.begin();
			for(std::vector< double>::iterator it = x_scores.begin(); it < x_scores.end(); ++it)
			{
				*probs = *it + minimum + 0.001;
				++probs;
			}
			sort(probabilities.begin(),probabilities.end());		
			DoubleReal negative_prior(0.5);
					
			GaussFitter::GaussFitResult	gauss_fit2,gauss_fit1;
			int x_score_start = ceil(probabilities.size()*0.7);
			gauss_fit1.x0 = gsl_stats_mean(&probabilities[x_score_start], 1, probabilities.size() -x_score_start) + probabilities[x_score_start];//(gauss_scores.begin()->getX() + (gauss_scores.end()-1)->getX())/2; 
			gauss_fit1.sigma = gsl_stats_sd(&probabilities[0], 1, probabilities.size() - 1);//pow(gsl_stats_sd_with_fixed_mean(&probabilities[x_score_start], 1, probabilities.size() - x_score_start, gauss_fit_param_.x0),2);
			gauss_fit1.A = 1	/sqrt(2*3.14159*pow(gauss_fit1.sigma,2));
					
			gauss_fit2.x0 = gsl_stats_mean(&probabilities[0], 1, ceil(0.5* probabilities.size())) + probabilities[0];
			gauss_fit2.sigma = gauss_fit1.sigma;
			gauss_fit2.A = 1	/sqrt(2*3.14159*pow(gauss_fit1.sigma,2));
			DoubleReal maxlike(0); 					
			for(vector<double >::const_iterator it = probabilities.begin(); it < probabilities.end() ; ++it)
			{
					DoubleReal the_x = *it;
					DoubleReal x_gauss2 = gauss_fit2.A * exp(-1.0 * pow(the_x - gauss_fit2.x0, 2) / (2 * pow(gauss_fit2.sigma, 2)));		
					DoubleReal x_gauss1 = gauss_fit1.A * exp(-1.0 * pow(the_x - gauss_fit1.x0, 2) / (2 * pow(gauss_fit1.sigma, 2)));
					maxlike += log10(negative_prior*x_gauss2+(1-negative_prior)*x_gauss1);
			}
			//-------------------------------------------------------------
			// create files for output
			//-------------------------------------------------------------
			bool output_plots  = param_.getValue("output_plots").toBool();
			String step;
			int iter = 0;
			String output;
			String output_ending;
			TextFile file;
			if(output_plots)
			{
				std::vector<DPosition<2> > points;
				DPosition<2> temp;
				double dividing_score = probabilities.back()/(Int)param_.getValue("number_of_bins");
				temp.setX(dividing_score/2);
				temp.setY(0);
				points.push_back(temp);
				double temp_divider = dividing_score;
				for(std::vector< double>::iterator it = probabilities.begin(); it < probabilities.end(); ++it)
				{
					if(temp_divider - *it >= 0)
					{
						points.back().setY(points.back().getY()+1);
					}
					else
					{
						temp.setX((temp_divider + temp_divider + dividing_score)/2);
						temp.setY(0);
						points.push_back(temp);
						temp_divider += dividing_score;
					}
				}	

				for(vector<DPosition<2> >::iterator it = points.begin(); it < points.end(); ++it)
				{
					it->setY( it->getY()/( probabilities.size() * dividing_score));
				}
					
				TextFile data_points;
				for(vector<DPosition<2> >::iterator it = points.begin(); it < points.end() ; ++it)
				{
					String temp  = it->getX();
					temp += "\t";
					temp += it->getY();
					data_points<<temp;
				}
				data_points.store((String)param_.getValue("output_name") + "scores.txt");
					
				step = (String)param_.getValue("output_name") + "step_";
				output = "set output \""+ (String)param_.getValue("output_name") +"step_";
				output_ending = ".pdf\"";
				iter = 0;
				file<<"set terminal pdf";	
				file<<(output + iter + output_ending);
				gauss_fit_param_ = gauss_fit2;
				String formula;
				formula =  "f(x)=" + String(gauss_fit_param_.A) +" * exp(-(x - " + String(gauss_fit_param_.x0) + ") ** 2 / 2 / (" + String(gauss_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior);
				file << formula;
				gauss_fit_param_ = gauss_fit1;
				formula = getGaussGnuplotFormula()+ "* (1 - " + String(negative_prior) + ")";
				file<<formula;
				formula = getBothGnuplotFormula(negative_prior); 
				file<<formula;
				file<<"plot \""+(String)param_.getValue("output_name") +"scores.txt\" with boxes, f(x) , g(x), h(x)";
				file.store(step + iter);	
			}
			//-------------------------------------------------------------
			// Estimate Parameters
			//-------------------------------------------------------------					
			bool stop_em_init = false;
			do
			{ 
				//E-STEP
				
				//sum new posterior probability
				DoubleReal sum_posterior(0), one_minus_sum_posterior(0);
				for(vector<double >::const_iterator it = probabilities.begin(); it < probabilities.end(); ++it)
				{
					DoubleReal the_x = *it;
					DoubleReal x_gauss2 = gauss_fit2.A * exp(-1.0 * pow(the_x - gauss_fit2.x0, 2) / (2 * pow(gauss_fit2.sigma, 2)));		
					DoubleReal x_gauss1 = gauss_fit1.A * exp(-1.0 * pow(the_x - gauss_fit1.x0, 2) / (2 * pow(gauss_fit1.sigma, 2)));
					sum_posterior += (negative_prior*x_gauss2)/((negative_prior*x_gauss2) + (1-negative_prior)*x_gauss1);
					one_minus_sum_posterior += 1  - ((negative_prior*x_gauss2)/((negative_prior*x_gauss2) + (1-negative_prior)*x_gauss1));
				}		
				//M-STEP
					
				//new mean
				DoubleReal sum_gauss1_x0(0),sum_gauss2_x0(0);
				for(vector<double >::const_iterator it = probabilities.begin(); it < probabilities.end() ; ++it)
				{
					DoubleReal the_x = *it;
					DoubleReal x_gauss2 = gauss_fit2.A * exp(-1.0 * pow(the_x - gauss_fit2.x0, 2) / (2 * pow(gauss_fit2.sigma, 2)));		
					DoubleReal x_gauss1 = gauss_fit1.A * exp(-1.0 * pow(the_x - gauss_fit1.x0, 2) / (2 * pow(gauss_fit1.sigma, 2)));
							
					sum_gauss1_x0 += ((1  - ((negative_prior*x_gauss2)/((negative_prior*x_gauss2) + (1-negative_prior)*x_gauss1)))*the_x);
					sum_gauss2_x0 += (((negative_prior*x_gauss2)/((negative_prior*x_gauss2) + (1-negative_prior)*x_gauss1))*the_x);
				}

										
				DoubleReal gauss1_mean = sum_gauss1_x0/one_minus_sum_posterior; 
				DoubleReal gauss2_mean = sum_gauss2_x0/sum_posterior; 
						
				//new standard deviation
				DoubleReal sum_gauss1_sigma(0),sum_gauss2_sigma(0);
				for(vector<double >::const_iterator it = probabilities.begin(); it < probabilities.end() ; ++it)
				{
					DoubleReal the_x = *it;
					DoubleReal x_gauss2 = gauss_fit2.A * exp(-1.0 * pow(the_x - gauss_fit2.x0, 2) / (2 * pow(gauss_fit2.sigma, 2)));		
					DoubleReal x_gauss1 = gauss_fit1.A * exp(-1.0 * pow(the_x - gauss_fit1.x0, 2) / (2 * pow(gauss_fit1.sigma, 2)));
							
					sum_gauss1_sigma += ((1  - ((negative_prior*x_gauss2)/((negative_prior*x_gauss2) + (1-negative_prior)*x_gauss1)))*pow(the_x - gauss1_mean,2));
					sum_gauss2_sigma += (((negative_prior*x_gauss2)/((negative_prior*x_gauss2) + (1-negative_prior)*x_gauss1))*pow(the_x - gauss2_mean,2));
				}					
						
				//update parameters
				gauss_fit1.x0 = gauss1_mean;
				gauss_fit1.sigma = sqrt(sum_gauss1_sigma/one_minus_sum_posterior); 
				gauss_fit1.A = 1/sqrt(2*3.14159*pow(gauss_fit1.sigma,2));
						
				gauss_fit2.x0 = gauss2_mean;
				gauss_fit2.sigma = sqrt(sum_gauss2_sigma/sum_posterior);
				gauss_fit2.A = 1/sqrt(2*3.14159*pow(gauss_fit2.sigma,2));

				//compute new prior probabilities negative peptides	
				sum_posterior = 0;
				for(vector<double >::const_iterator it = probabilities.begin(); it < probabilities.end() ; ++it)
				{
					DoubleReal the_x = *it;
					DoubleReal x_gauss2 = gauss_fit2.A * exp(-1.0 * pow(the_x - gauss_fit2.x0, 2) / (2 * pow(gauss_fit2.sigma, 2)));		
					DoubleReal x_gauss1 = gauss_fit1.A * exp(-1.0 * pow(the_x - gauss_fit1.x0, 2) / (2 * pow(gauss_fit1.sigma, 2)));
					sum_posterior += (negative_prior*x_gauss2)/((negative_prior*x_gauss2) + (1-negative_prior)*x_gauss1);
				}
						
				negative_prior = sum_posterior/probabilities.size();
						
           	
        //compute new maximum likelihood
        DoubleReal new_maxlike(0);
        for(vector<double >::const_iterator it = probabilities.begin(); it < probabilities.end() ; ++it)
				{
					DoubleReal the_x = *it;
					DoubleReal x_gauss2 = gauss_fit2.A * exp(-1.0 * pow(the_x - gauss_fit2.x0, 2) / (2 * pow(gauss_fit2.sigma, 2)));		
					DoubleReal x_gauss1 = gauss_fit1.A * exp(-1.0 * pow(the_x - gauss_fit1.x0, 2) / (2 * pow(gauss_fit1.sigma, 2)));
					new_maxlike += log10(negative_prior*x_gauss2+(1-negative_prior)*x_gauss1);
				}
        if(fabs(new_maxlike - maxlike) < 0.001)
        {
        	stop_em_init = true;      		
        }
				if(output_plots)
				{
					++iter;
					file[1] = (output + iter + output_ending);
					gauss_fit_param_ = gauss_fit2;
					String formula;
					formula =  "f(x)=" + String(gauss_fit_param_.A) + " * exp(-(x - " + String(gauss_fit_param_.x0) + ") ** 2 / 2 / (" + String(gauss_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior);
					file[2] = formula;
					gauss_fit_param_ = gauss_fit1;
					file[3] = getGaussGnuplotFormula()+ "* (1 - " + String(negative_prior) + ")"; 
					formula = "h(x)=" + String(negative_prior)+"*" + String(gauss_fit2.A) + " * exp(-(x - " + String(gauss_fit2.x0) + ") ** 2 / 2 / (" + String(gauss_fit2.sigma) + ") ** 2)"+ " + "+"(1-"+String(negative_prior) + ")*" + String(gauss_fit1.A) + " * exp(-(x - " + String(gauss_fit1.x0) + ") ** 2 / 2 / (" + String(gauss_fit1.sigma) + ") ** 2)";
					file[4] = formula;
					file.store(step + iter);
				}
				//update maximum likelihood
				maxlike = new_maxlike;
			}while(!stop_em_init);									
					
			if(gauss_fit1.x0 <= gauss_fit2.x0)
			{
				gauss_fit_param_ = gauss_fit2;
				gumbel_fit_param_.b =gauss_fit1.sigma;// 6*gauss_fit1.sigma/pow(3.14,2);
				gumbel_fit_param_.a = gauss_fit1.x0; // gauss_fit1.x0 + 0.57722* 1/(3.14*sqrt(6*gumbel_fit_param_.b));                         
			}
			else
			{
				gauss_fit_param_ = gauss_fit1;
				gumbel_fit_param_.b =gauss_fit2.sigma; //6*gauss_fit2.sigma/pow(3.14,2);
				gumbel_fit_param_.a =gauss_fit2.x0; //gauss_fit2.x0 + 0.57722* 1/(3.14*sqrt(6*gumbel_fit_param_.b));
			}
			if(output_plots)
			{
				++iter;			
				file[1] = (output + iter + output_ending);
				file[2] = getGumbelGnuplotFormula() + "*" + String(negative_prior);
				file[3] = getGaussGnuplotFormula()+ "* (1 - " + String(negative_prior) + ")";
				file[4] = getBothGnuplotFormula(negative_prior);
				file.store(step + iter);				
			}
			//Compute probabilities
			max_gumbel_ = exp(-1.0)/gumbel_fit_param_.b;
			max_gauss_ = gauss_fit_param_.A;
			negative_prior_ = negative_prior;
			smallest_score_ = minimum;
			probs = probabilities.begin();
			for(vector<double >::iterator it = x_scores.begin(); it < x_scores.end() ; ++it, ++probs)
			{
				*probs = 	computeProbability(*it);			
			}
		}
				
		const String PosteriorErrorProbabilityModel::getGumbelGnuplotFormula() const
		{
			// build a formula with the fitted parameters for gnuplot
			stringstream formula;
			formula <<"f(x)=" << "(1/" << gumbel_fit_param_.b <<") * " << "exp(( "<< gumbel_fit_param_.a<< "- x)/"<< gumbel_fit_param_.b <<") * exp(-exp(("<<gumbel_fit_param_.a<<" - x)/"<<gumbel_fit_param_.b<<"))";
			return  formula.str();
		}
				
		const String PosteriorErrorProbabilityModel::getGaussGnuplotFormula() const
		{
			stringstream formula;
			formula << "g(x)=" << gauss_fit_param_.A << " * exp(-(x - " << gauss_fit_param_.x0 << ") ** 2 / 2 / (" << gauss_fit_param_.sigma << ") ** 2)";
			return formula.str();
		}
				
		const String PosteriorErrorProbabilityModel::getBothGnuplotFormula(double negative_prior) const
		{
			stringstream formula;
			formula << "h(x)=" << negative_prior<<"*"<< "(1/" << gumbel_fit_param_.b <<") * " << "exp(( "<< gumbel_fit_param_.a<< "- x)/"<< gumbel_fit_param_.b <<") * exp(-exp(("<<gumbel_fit_param_.a<<" - x)/"<<gumbel_fit_param_.b<<")) + "<<"(1-"<<negative_prior<<")*"<<gauss_fit_param_.A << " * exp(-(x - " << gauss_fit_param_.x0 << ") ** 2 / 2 / (" << gauss_fit_param_.sigma << ") ** 2)";
			return formula.str();
		}	
		
		DoubleReal PosteriorErrorProbabilityModel::computeProbability(DoubleReal score)
		{
			DoubleReal the_x = score + smallest_score_ + 0.001;
			DoubleReal x_gumbel;
			DoubleReal x_gauss;
			if(the_x < gumbel_fit_param_.a)
			{
				x_gumbel = max_gumbel_;	
				x_gauss = gauss_fit_param_.A * exp(-1.0 * pow(the_x - gauss_fit_param_.x0, 2) / (2 * pow(gauss_fit_param_.sigma, 2)));
			}
			else if(the_x > gauss_fit_param_.x0)
			{
				DoubleReal z = exp((gumbel_fit_param_.a - the_x)/gumbel_fit_param_.b);
				x_gumbel = (z*exp(-1* z))/gumbel_fit_param_.b;				
				x_gauss = max_gauss_;
			}
			else
			{
				DoubleReal z = exp((gumbel_fit_param_.a - the_x)/gumbel_fit_param_.b);
				x_gumbel = (z*exp(-1* z))/gumbel_fit_param_.b;				
				x_gauss = gauss_fit_param_.A * exp(-1.0 * pow(the_x - gauss_fit_param_.x0, 2) / (2 * pow(gauss_fit_param_.sigma, 2)));
			}
			 return (negative_prior_*x_gumbel)/((negative_prior_*x_gumbel) + (1-negative_prior_)*x_gauss);		
		}
		
	} //namespace Math
} // namespace OpenMS

