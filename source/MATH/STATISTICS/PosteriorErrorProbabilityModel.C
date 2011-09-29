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
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
//
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <algorithm>
#include <gsl/gsl_statistics.h>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

namespace OpenMS
{
	namespace Math
	{
		PosteriorErrorProbabilityModel::PosteriorErrorProbabilityModel()
		: DefaultParamHandler("PosteriorErrorProbabilityModel"), negative_prior_(0.5),max_incorrectly_(0),max_correctly_(0), smallest_score_(0)
  	{
			defaults_.setValue("number_of_bins", 100, "Number of bins used for visualization. Only needed if each iteration step of the EM-Algorithm will be visualized", StringList::create("advanced"));
			defaults_.setValue("output_plots","false","If true every step of the EM-algorithm will be written to a file as a gnuplot formula",StringList::create("advanced"));
			defaults_.setValidStrings("output_plots",StringList::create("true,false"));
			defaults_.setValue("output_name","", "if output_plots is on, the output files will be saved in the following manner: <output_name>scores.txt for the scores and <output_name> which contains each step of the EM-algorithm e.g. output_name = /usr/home/OMSSA123 then /usr/home/OMSSA123_scores.txt, /usr/home/OMSSA123 will be written. If no directory is specified, e.g. instead of '/usr/home/OMSSA123' just OMSSA123, the files will be written into the working directory.",StringList::create("advanced,output file"));
			defaults_.setValue("incorrectly_assigned","Gumbel", "for 'Gumbel', the Gumbel distribution is used to plot incorrectly assigned sequences. For 'Gauss', the Gauss distribution is used.",StringList::create("advanced"));
			defaults_.setValidStrings("incorrectly_assigned",StringList::create("Gumbel,Gauss"));
			defaultsToParam_();
			calc_incorrect_ = &PosteriorErrorProbabilityModel::getGumbel;
			calc_correct_ = &PosteriorErrorProbabilityModel::getGauss;
			getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGumbelGnuplotFormula;
			getPositiveGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;
		}
		
		PosteriorErrorProbabilityModel::~PosteriorErrorProbabilityModel()
		{
		}
				
		void PosteriorErrorProbabilityModel::fit( std::vector<double>& search_engine_scores)
		{	
			if(search_engine_scores.empty())
			{
				return;
			}
			//-------------------------------------------------------------
			// Initializing Parameters
			//-------------------------------------------------------------
			sort(search_engine_scores.begin(),search_engine_scores.end());		
			
			smallest_score_ = search_engine_scores[0];
			vector<double> x_scores;
			x_scores.resize(search_engine_scores.size());
			std::vector< double>::iterator it = x_scores.begin();
			for(std::vector< double>::iterator iti = search_engine_scores.begin(); iti < search_engine_scores.end(); ++it,++iti)
			{
				*it = *iti + fabs(smallest_score_) + 0.001;
			}
			negative_prior_ = 0.7;
			if(param_.getValue("incorrectly_assigned") == "Gumbel")
			{
				incorrectly_assigned_fit_param_.x0 = gsl_stats_mean(&x_scores[0], 1, ceil(0.5* x_scores.size())) + x_scores[0];
				incorrectly_assigned_fit_param_.sigma = gsl_stats_sd(&x_scores[0], 1, x_scores.size() - 1);//pow(gsl_stats_sd_with_fixed_mean(&probabilities[x_score_start], 1, probabilities.size() - x_score_start, gauss_fit_param_.x0),2);
				incorrectly_assigned_fit_param_.A = 1	/sqrt(2*3.14159*pow(incorrectly_assigned_fit_param_.sigma,2));	
				//TODO: compute directly with gauss. Workaround:
				calc_incorrect_ = &PosteriorErrorProbabilityModel::getGauss;
				getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGumbelGnuplotFormula;
			}
			else
			{
				incorrectly_assigned_fit_param_.x0 = gsl_stats_mean(&x_scores[0], 1, ceil(0.5* x_scores.size())) + x_scores[0];
				incorrectly_assigned_fit_param_.sigma = gsl_stats_sd(&x_scores[0], 1, x_scores.size() - 1);//pow(gsl_stats_sd_with_fixed_mean(&probabilities[x_score_start], 1, probabilities.size() - x_score_start, gauss_fit_param_.x0),2);
				incorrectly_assigned_fit_param_.A = 1	/sqrt(2*3.14159*pow(incorrectly_assigned_fit_param_.sigma,2));		
				calc_incorrect_ = &PosteriorErrorProbabilityModel::getGauss;
				getNegativeGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;
			}
			getPositiveGnuplotFormula_ = &PosteriorErrorProbabilityModel::getGaussGnuplotFormula;
			calc_correct_ = &PosteriorErrorProbabilityModel::getGauss;
			Size x_score_start = std::min(x_scores.size()-1, (Size) ceil(x_scores.size()*0.7)); // if only one score is present, ceil(...) will yield 1, which is an invalid index
			correctly_assigned_fit_param_.x0 = gsl_stats_mean(&x_scores[x_score_start], 1, x_scores.size() - x_score_start) + x_scores[x_score_start];//(gauss_scores.begin()->getX() + (gauss_scores.end()-1)->getX())/2; 
			correctly_assigned_fit_param_.sigma = incorrectly_assigned_fit_param_.sigma;
			correctly_assigned_fit_param_.A = 1.0	/ sqrt(2*3.14159*pow(correctly_assigned_fit_param_.sigma,2));		

			DoubleReal maxlike(0);
			vector<DoubleReal> incorrect_density;
			vector<DoubleReal> correct_density;
			
			fillDensities(x_scores,incorrect_density,correct_density);
			
			
			maxlike = computeMaxLikelihood(incorrect_density,correct_density);
			//-------------------------------------------------------------
			// create files for output
			//-------------------------------------------------------------
			bool output_plots  = param_.getValue("output_plots").toBool();
			TextFile* file = NULL;
			if(output_plots)
			{
				file = InitPlots(x_scores);
			}			
			//-------------------------------------------------------------
			// Estimate Parameters - EM algorithm
			//-------------------------------------------------------------
			bool stop_em_init = false;
			do
			{ 
				//E-STEP
				DoubleReal one_minus_sum_posterior = one_minus_sum_post(incorrect_density,correct_density);
				DoubleReal sum_posterior = sum_post(incorrect_density,correct_density);
				
				//new mean				
				DoubleReal sum_positive_x0 = sum_pos_x0(x_scores, incorrect_density,correct_density);
				DoubleReal sum_negative_x0 = sum_neg_x0(x_scores, incorrect_density,correct_density);
				
				DoubleReal positive_mean = sum_positive_x0/one_minus_sum_posterior; 
				DoubleReal negative_mean = sum_negative_x0/sum_posterior; 
				
				//new standard deviation
				DoubleReal sum_positive_sigma = sum_pos_sigma(x_scores,incorrect_density,correct_density, positive_mean);
				DoubleReal sum_negative_sigma = sum_neg_sigma(x_scores,incorrect_density,correct_density, negative_mean);
				
				//update parameters
				correctly_assigned_fit_param_.x0 = positive_mean;
				if(sum_positive_sigma  != 0 && one_minus_sum_posterior != 0)
				{
					correctly_assigned_fit_param_.sigma = sqrt(sum_positive_sigma/one_minus_sum_posterior); 
					correctly_assigned_fit_param_.A = 1/sqrt(2*3.14159*pow(correctly_assigned_fit_param_.sigma,2));
				}
						
				incorrectly_assigned_fit_param_.x0 = negative_mean;
				if(sum_negative_sigma  != 0 && sum_posterior != 0)
				{
					incorrectly_assigned_fit_param_.sigma = sqrt(sum_negative_sigma/sum_posterior);
					incorrectly_assigned_fit_param_.A = 1/sqrt(2*3.14159*pow(incorrectly_assigned_fit_param_.sigma,2));
				}				

					
				//compute new prior probabilities negative peptides
				fillDensities(x_scores,incorrect_density,correct_density);
				sum_posterior = sum_post(incorrect_density,correct_density);
				negative_prior_ = sum_posterior/x_scores.size();
				
				DoubleReal new_maxlike(computeMaxLikelihood(incorrect_density,correct_density));
        if(boost::math::isnan(new_maxlike - maxlike))
				{
					throw Exception::UnableToFit(__FILE__,__LINE__,__PRETTY_FUNCTION__,"UnableToFit-PosteriorErrorProbability","Could not fit mixture model to data");					
				}
        if(fabs(new_maxlike - maxlike) < 0.001)
        {
        	stop_em_init = true;
        	sum_posterior = sum_post(incorrect_density,correct_density);
					negative_prior_ = sum_posterior/x_scores.size();
        	
        }
				if(output_plots)
				{
					String formula1, formula2,formula3;
					formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_)+ "* " + String(negative_prior_);//String(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + String(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior_);
					formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_)+ "* (1 - " + String(negative_prior_) + ")";//String(correctly_assigned_fit_param_.A) +" * exp(-(x - " + String(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + String(negative_prior_) + ")";
					formula3 = getBothGnuplotFormula(incorrectly_assigned_fit_param_,correctly_assigned_fit_param_);
					(*file)<<("plot \""+(String)param_.getValue("output_name") +"_scores.txt\" with boxes, " + formula1 + " , " + formula2 + " , " + formula3);
				}
				//update maximum likelihood
				maxlike = new_maxlike;				
			}while(!stop_em_init);
			//-------------------------------------------------------------
			// Finished fitting
			//-------------------------------------------------------------			
			//!!Workaround:
			if(param_.getValue("incorrectly_assigned") == "Gumbel")
			{
				calc_incorrect_ = &PosteriorErrorProbabilityModel::getGumbel;
			}
			max_incorrectly_ = ((this)->*(calc_incorrect_))(incorrectly_assigned_fit_param_.x0, incorrectly_assigned_fit_param_);
			max_correctly_ = ((this)->*(calc_correct_))(correctly_assigned_fit_param_.x0, correctly_assigned_fit_param_ );
			if(output_plots)
			{
				String formula1, formula2,formula3;
				formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_)+ "*" + String(negative_prior_);//String(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + String(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior_);
				formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_)+ "* (1 - " + String(negative_prior_) + ")";// String(correctly_assigned_fit_param_.A) +" * exp(-(x - " + String(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + String(negative_prior_) + ")";
				formula3 = getBothGnuplotFormula(incorrectly_assigned_fit_param_,correctly_assigned_fit_param_);
				(*file)<<("plot \""+(String)param_.getValue("output_name") +"_scores.txt\" with boxes, " + formula1 + " , " + formula2 + " , " + formula3);
				file->store((String)param_.getValue("output_name"));
				delete file;
			}
		}
		
		void PosteriorErrorProbabilityModel::fit(  std::vector<double>& search_engine_scores, vector<double>& probabilities)
		{	
			fit(search_engine_scores);
			probabilities.resize(search_engine_scores.size());
			vector<double>::iterator probs =probabilities.begin();
			for(vector<double>::iterator scores = search_engine_scores.begin(); scores != search_engine_scores.end(); ++scores, ++probs)
			{
				*probs = computeProbability(*scores);
			}
		}
		
		void PosteriorErrorProbabilityModel::fillDensities(vector<double>& x_scores,vector<DoubleReal>& incorrect_density,vector<DoubleReal>& correct_density)
		{
			if(incorrect_density.size() != x_scores.size())
			{
				incorrect_density.resize(x_scores.size());
				correct_density.resize(x_scores.size());
			}
			vector<DoubleReal>::iterator incorrect = incorrect_density.begin();
			vector<DoubleReal>::iterator correct = correct_density.begin();
			for(vector<double>::iterator scores = x_scores.begin(); scores != x_scores.end(); ++scores, ++incorrect, ++correct)
			{
				*incorrect = ((this)->*(calc_incorrect_))(*scores, incorrectly_assigned_fit_param_);
				*correct = ((this)->*(calc_correct_))(*scores, correctly_assigned_fit_param_ );
			}
		}
		
		DoubleReal PosteriorErrorProbabilityModel::computeMaxLikelihood(vector<DoubleReal>& incorrect_density, vector<DoubleReal>& correct_density)
		{
			DoubleReal maxlike(0);
			vector<DoubleReal>::iterator incorrect = incorrect_density.begin();
			for(vector<DoubleReal>::iterator correct = correct_density.begin(); correct < correct_density.end() ; ++correct, ++incorrect)
			{
				maxlike += log10(negative_prior_* (*incorrect) +(1-negative_prior_)* (*correct));
			}
			return maxlike;
		}
		DoubleReal PosteriorErrorProbabilityModel::one_minus_sum_post(vector<DoubleReal>& incorrect_density, vector<DoubleReal>& correct_density)
		{
			DoubleReal one_min(0);
			vector<DoubleReal>::iterator incorrect = incorrect_density.begin();
			for(vector<DoubleReal>::iterator correct = correct_density.begin(); correct < correct_density.end() ; ++correct, ++incorrect)
			{
				one_min +=  1  - ((negative_prior_* (*incorrect) )/((negative_prior_* (*incorrect) ) + (1-negative_prior_)* (*correct) ));
			}
			return one_min;		
		}

		DoubleReal PosteriorErrorProbabilityModel::sum_post(vector<DoubleReal>& incorrect_density, vector<DoubleReal>& correct_density)
		{
			DoubleReal post(0);
			vector<DoubleReal>::iterator incorrect = incorrect_density.begin();
			for(vector<DoubleReal>::iterator correct = correct_density.begin(); correct < correct_density.end() ; ++correct, ++incorrect)
			{
				post += ((negative_prior_* (*incorrect) )/((negative_prior_* (*incorrect) ) + (1-negative_prior_)* (*correct) ));
			}
			return post;		
		}
		
		DoubleReal PosteriorErrorProbabilityModel::sum_pos_x0(vector<double>& x_scores, vector<DoubleReal>& incorrect_density, vector<DoubleReal>& correct_density)
		{
			DoubleReal pos_x0(0);
			vector<double>::iterator the_x = x_scores.begin();
			vector<DoubleReal>::iterator incorrect = incorrect_density.begin();
			for(vector<DoubleReal>::iterator correct = correct_density.begin(); correct < correct_density.end() ; ++correct, ++incorrect,++the_x)
			{
				pos_x0 += ((1  - ((negative_prior_* (*incorrect) )/((negative_prior_* (*incorrect) ) + (1-negative_prior_)* (*correct) )))* (*the_x) );
			}
			return pos_x0;	
		}
		
		DoubleReal PosteriorErrorProbabilityModel::sum_neg_x0(vector<double>& x_scores, vector<DoubleReal>& incorrect_density, vector<DoubleReal>& correct_density)
		{
			DoubleReal neg_x0(0);
			vector<double>::iterator the_x = x_scores.begin();
			vector<DoubleReal>::iterator correct = correct_density.begin();
			for(vector<DoubleReal>::iterator incorrect = incorrect_density.begin(); incorrect < incorrect_density.end() ; ++correct, ++incorrect,++the_x)
			{
				neg_x0 += ((((negative_prior_* (*incorrect) )/((negative_prior_* (*incorrect) ) + (1-negative_prior_)* (*correct) )))* (*the_x) );
			}
			return neg_x0;	
		}		
		
		DoubleReal PosteriorErrorProbabilityModel::sum_pos_sigma(vector<double>& x_scores, vector<DoubleReal>& incorrect_density, vector<DoubleReal>& correct_density, DoubleReal positive_mean)
		{
			DoubleReal pos_sigma(0);
			vector<double>::iterator the_x = x_scores.begin();
			vector<DoubleReal>::iterator incorrect = incorrect_density.begin();
			for(vector<DoubleReal>::iterator correct = correct_density.begin(); correct < correct_density.end() ; ++correct, ++incorrect,++the_x)
			{
				pos_sigma += ((1  - ((negative_prior_* (*incorrect) )/((negative_prior_* (*incorrect) ) + (1-negative_prior_)* (*correct) )))*pow( (*the_x) - positive_mean,2));
			}
			return pos_sigma;
		}
		
		DoubleReal PosteriorErrorProbabilityModel::sum_neg_sigma(vector<double>& x_scores, vector<DoubleReal>& incorrect_density, vector<DoubleReal>& correct_density, DoubleReal positive_mean)
		{
			DoubleReal neg_sigma(0);
			vector<double>::iterator the_x = x_scores.begin();
			vector<DoubleReal>::iterator incorrect = incorrect_density.begin();
			for(vector<DoubleReal>::iterator correct = correct_density.begin(); correct < correct_density.end() ; ++correct, ++incorrect,++the_x)
			{
				neg_sigma += ((((negative_prior_* (*incorrect) )/((negative_prior_* (*incorrect) ) + (1-negative_prior_)* (*correct) )))*pow( (*the_x) - positive_mean,2));
			}
			return neg_sigma;
		}
		DoubleReal PosteriorErrorProbabilityModel::computeProbability(DoubleReal score)
		{
			score = score + fabs(smallest_score_) + 0.001;
			DoubleReal x_neg;
			DoubleReal x_pos;
			//the score is smaller than the peak of incorreclty assigned sequences. To ensure that the probabilies wont rise again use the incorrectly assigend peak for computation
			if(score < incorrectly_assigned_fit_param_.x0)
			{
				x_neg = max_incorrectly_;	
				x_pos = ((this)->*(calc_correct_))(score, correctly_assigned_fit_param_ );
			}
			//same as above. However, this time to ensure that probabilities wont drop again.
			else if(score > correctly_assigned_fit_param_.x0)
			{
				x_neg = ((this)->*(calc_incorrect_))(score, incorrectly_assigned_fit_param_);
				x_pos = max_correctly_;
			}
			//if its in between use the normal formula
			else
			{
				x_neg = ((this)->*(calc_incorrect_))(score, incorrectly_assigned_fit_param_);
				x_pos = ((this)->*(calc_correct_))(score, correctly_assigned_fit_param_ );
			}
			 return (negative_prior_*x_neg)/((negative_prior_*x_neg) + (1-negative_prior_)*x_pos);		
		}		
		TextFile* PosteriorErrorProbabilityModel::InitPlots(vector<double> & x_scores)
		{		
				TextFile* file = new TextFile;
				String output;
				std::vector<DPosition<2> > points;
				Int number_of_bins = param_.getValue("number_of_bins");
				points.resize(number_of_bins);
				DPosition<2> temp;
				double dividing_score = (x_scores.back() -x_scores[0])/number_of_bins;

				temp.setX(dividing_score/2);
				temp.setY(0);
				Int bin = 0;
				points[bin] = temp;
				double temp_divider = dividing_score;
				for(std::vector< double>::iterator it = x_scores.begin(); it < x_scores.end(); ++it)
				{
					if(temp_divider - *it >= 0 && bin < number_of_bins - 1 )
					{
							points[bin].setY(points[bin].getY()+1);
					}
					else if(bin  == number_of_bins - 1 )
					{
						points[bin].setY(points[bin].getY()+1);
					}
					else
					{
						temp.setX((temp_divider + temp_divider + dividing_score)/2);
						temp.setY(1);
						++bin;
						points[bin] = temp;
						temp_divider += dividing_score;
					}
				}	

				for(vector<DPosition<2> >::iterator it = points.begin(); it < points.end(); ++it)
				{
					it->setY( it->getY()/( x_scores.size()  *dividing_score));
				}
					
				TextFile data_points;
				for(vector<DPosition<2> >::iterator it = points.begin(); it < points.end() ; ++it)
				{
					String temp  = it->getX();
					temp += "\t";
					temp += it->getY();
					data_points<<temp;
				}
				data_points.store((String)param_.getValue("output_name") + "_scores.txt");
				output = "set output \""+ (String)param_.getValue("output_name") +".ps\"";
				(*file)<<"set terminal postscript color solid linewidth 2.0 rounded";
				//(*file)<<"set style empty solid 0.5 border -1";
				//(*file)<<"set style function lines";
				(*file)<<"set xlabel \"discriminant score\"";
				(*file)<<"set ylabel \"density\"";
				//TODO: (*file)<<"set title ";
				(*file)<<"set key off";
				(*file)<<(output);
				String formula1, formula2;
				formula1 = ((this)->*(getNegativeGnuplotFormula_))(incorrectly_assigned_fit_param_)+ "* " + String(negative_prior_);//String(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + String(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior_);
				formula2 = ((this)->*(getPositiveGnuplotFormula_))(correctly_assigned_fit_param_)+ "* (1 - " + String(negative_prior_) + ")";//String(correctly_assigned_fit_param_.A) +" * exp(-(x - " + String(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + String(negative_prior_) + ")";
				(*file)<< ("plot \""+(String)param_.getValue("output_name") +"_scores.txt\" with boxes, " + formula1 + " , " + formula2);
				return file;
		}

		const String PosteriorErrorProbabilityModel::getGumbelGnuplotFormula(const GaussFitter::GaussFitResult& params) const
		{
			// build a formula with the fitted parameters for gnuplot
			stringstream formula;
			formula << "(1/" << params.sigma <<") * " << "exp(( "<< params.x0<< "- x)/"<< params.sigma <<") * exp(-exp(("<<params.x0<<" - x)/"<<params.sigma<<"))";
			return  formula.str();
		}
				
		const String PosteriorErrorProbabilityModel::getGaussGnuplotFormula(const GaussFitter::GaussFitResult& params) const
		{
			stringstream formula;
			formula << params.A << " * exp(-(x - " << params.x0 << ") ** 2 / 2 / (" << params.sigma << ") ** 2)";
			return formula.str();
		}
		
		const String PosteriorErrorProbabilityModel::getBothGnuplotFormula(const GaussFitter::GaussFitResult& incorrect,const GaussFitter::GaussFitResult& correct) const
		{
			stringstream formula;
			formula << negative_prior_<<"*"<<  ((this)->*(getNegativeGnuplotFormula_))(incorrect) <<" + (1-"<<negative_prior_<<")*"<< ((this)->*(getPositiveGnuplotFormula_))(correct);
			return formula.str();
		}
		
		void PosteriorErrorProbabilityModel::plotTargetDecoyEstimation(vector<double> &target,vector<double> & decoy)
		{
				TextFile file;
				String output;
				std::vector<DPosition<3> > points;
				Int number_of_bins = param_.getValue("number_of_bins");
				points.resize(number_of_bins);
				DPosition<3> temp;

				sort(target.begin(),target.end());
				sort(decoy.begin(),decoy.end());
				
				double dividing_score = (max(target.back(),decoy.back())/*scores.back()*/ - min(target[0],decoy[0])/*scores[0]*/)/number_of_bins;
				
				temp[0] = (dividing_score/2);
				temp[1] = 0;
				temp[2] = 0;
				Int bin = 0;
				points[bin] = temp;
				double temp_divider = dividing_score;
				for(std::vector< double>::iterator it = target.begin(); it < target.end(); ++it)
				{
					*it = *it + fabs(smallest_score_) + 0.001;
					if(temp_divider - *it >= 0 && bin < number_of_bins - 1 )
					{
							points[bin][1] = (points[bin][1]+1);
					}
					else if(bin  == number_of_bins - 1 )
					{
						points[bin][1] = (points[bin][1]+1);
					}
					else
					{
						temp[0] = ((temp_divider + temp_divider + dividing_score)/2);
						temp[1] = 1;
						++bin;
						points[bin] = temp;
						temp_divider += dividing_score;
					}
				}	

				bin = 0;
				temp_divider = dividing_score;
				for(std::vector< double>::iterator it = decoy.begin(); it < decoy.end(); ++it)
				{	
					*it = *it + fabs(smallest_score_) + 0.001;
					if(temp_divider - *it >= 0 && bin < number_of_bins - 1 )
					{
							points[bin][2] = (points[bin][2]+1);
					}
					else if(bin  == number_of_bins - 1 )
					{
						points[bin][2] = (points[bin][2]+1);
					}
					else
					{
						//temp[0] = ((temp_divider + temp_divider + dividing_score)/2);
					//	temp[2] = 1;
						++bin;
						points[bin][2] = 1;
						temp_divider += dividing_score;
					}
				}	

				for(vector<DPosition<3> >::iterator it = points.begin(); it < points.end(); ++it)
				{
				//	if((*it)[1] > (*it)[2])
				//	{(*it)[1] = (*it)[1] + (*it)[2];}
				/*	else{/(*it)[2] = (*it)[1] + (*it)[2];//}*/
					
					(*it)[1] = ( (*it)[1]/( (decoy.size()+target.size())  *dividing_score));
					(*it)[2] = ( (*it)[2]/( (decoy.size()+target.size())  *dividing_score));
				}
					
				TextFile data_points;
				for(vector<DPosition<3> >::iterator it = points.begin(); it < points.end() ; ++it)
				{
					String temp  = (*it)[0];
					temp += "\t";
					temp += (*it)[1];
					temp += "\t";
					temp += (*it)[2];
					data_points<<temp;
				}
				data_points.store( (String)param_.getValue("output_name") + "_target_decoy_scores.txt");
				output = String("set output \"")+  (String)param_.getValue("output_name") + "_target_decoy.ps\"";
				(file)<<"set terminal postscript color solid linewidth 2.0 rounded";
				//(*file)<<"set style empty solid 0.5 border -1";
				//(*file)<<"set style function lines";
				(file)<<"set xlabel \"discriminant score\"";
				(file)<<"set ylabel \"density\"";
				//TODO: (*file)<<"set title ";
				(file)<<"set key off";
				(file)<<(output);
				String formula1, formula2;
				formula1 = getGumbelGnuplotFormula(getIncorrectlyAssignedFitResult())+ "* " + String(getNegativePrior());//String(incorrectly_assigned_fit_param_.A) +" * exp(-(x - " + String(incorrectly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(incorrectly_assigned_fit_param_.sigma) + ") ** 2)"+ "*" + String(negative_prior_);
				formula2 = getGaussGnuplotFormula(getCorrectlyAssignedFitResult())+ "* (1 - " + String(getNegativePrior()) + ")";//String(correctly_assigned_fit_param_.A) +" * exp(-(x - " + String(correctly_assigned_fit_param_.x0) + ") ** 2 / 2 / (" + String(correctly_assigned_fit_param_.sigma) + ") ** 2)"+ "* (1 - " + String(negative_prior_) + ")";
				(file)<< ("plot \"" + (String)param_.getValue("output_name") + "_target_decoy_scores.txt\"   using 1:3  with boxes fill solid 0.8 noborder, \"" + (String)param_.getValue("output_name") + "_target_decoy_scores.txt\"  using 1:2  with boxes, " + formula1 + " , " + formula2);
				file.store( (String)param_.getValue("output_name") + "_target_decoy");
			}		
		
	} //namespace Math
} // namespace OpenMS

