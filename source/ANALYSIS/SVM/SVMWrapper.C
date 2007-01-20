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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <numeric>
#include <iostream>
#include <fstream>
#include <ctime>

#include <gsl/gsl_cdf.h>

using namespace std;

namespace OpenMS 
{
		
	SVMWrapper::SVMWrapper() : param_(NULL), model_(NULL), sigma_(0), sigmas_(vector<DoubleReal>()), gauss_table_(), kernel_type_(PRECOMPUTED), border_length_(0), training_set_(NULL), training_problem_(NULL)
	{	
	  param_ = (struct svm_parameter*) malloc(sizeof(struct svm_parameter));
	  initParameters();
	}
	
	SVMWrapper::~SVMWrapper()
	{
	  if (param_ != NULL)
	  {
			free(param_);
			param_ = NULL;
	  }
	  if (model_ != NULL)
	  {
			svm_destroy_model(model_);
			model_ = NULL;
	  } 
	}
	
	void SVMWrapper::setParameter(SVM_parameter_type type, int value)
	{
	
	  switch(type)
	  {
	  	case(SVM_TYPE):
				if (value == NU_SVR || value == EPSILON_SVR || value == NU_SVC || value == C_SVC || value == ONE_CLASS)	
				{
				  param_->svm_type = value;
				}
				break;
			case(KERNEL_TYPE):
				kernel_type_ = value;
				if (value == OLIGO)
				{
					param_->kernel_type = PRECOMPUTED;
				}
				else
				{	
					param_->kernel_type = value;
				}
				break;
			case(DEGREE):
				param_->degree = value;
			  break;
			case(BORDER_LENGTH):
				border_length_ = value;
			  break;
			case(C):
				param_->C = value;
				break;
			case(P):
				param_->p = value;
				break;
			case(NU):
				param_->nu = value;
				break;
			case(GAMMA):
				param_->gamma = value;
				break;
			case(SIGMA):
				sigma_ = value;
				if (border_length_ >= 1)
				{
					SVMWrapper::calculateGaussTable(border_length_, sigma_, gauss_table_);
    		}
				
				break;
			case(PROBABILITY):
				if (value == 1 || value == 0)
				{
					param_->probability = value;
				}
				break;
			default:
			  break;
	  }
	}
	
	int SVMWrapper::getIntParameter(SVM_parameter_type type)
	{
	
	    switch(type)
	    {
	    	case(KERNEL_TYPE):
					if (param_->kernel_type != PRECOMPUTED)
					{
	    	    return param_->kernel_type;
	    	    break;
					}
					else
					{
						return kernel_type_;
						break;
					}
	    	case(SVM_TYPE):
	    	    return param_->svm_type;
	    	    break;
	    	case(DEGREE):
	    	    return (int) param_->degree;
	    	    break;
	    	case(BORDER_LENGTH):
	    	    return border_length_;
	    	    break;
	    	case(PROBABILITY):
	    	    return param_->probability;
	    	    break;
	    	default:
	          return -1;
	    	    break;
	    }
	}
	
	void SVMWrapper::setParameter(SVM_parameter_type type, double value)
	{
	    switch(type)
	    {
				case(DEGREE):
					param_->degree = (int)value;
				  break;
				case(C):
			    param_->C = value;
		   		break;
				case(P):
			    param_->p = value;
		   		break;
				case(NU):
			    param_->nu = value;
		   		break;
				case(GAMMA):
			    param_->gamma = value;
		   		break;
				case(SIGMA):
					sigma_ = value;
					if (border_length_ >= 1)
					{
						SVMWrapper::calculateGaussTable(border_length_, sigma_, gauss_table_);
	    		}
					break;
				default:
			    break;
	    }
	}
	
	double SVMWrapper::getDoubleParameter(SVM_parameter_type type)
	{
    switch(type)
    {
     	case(C):
     	    return param_->C;
     	    break;
     	case(P):
     	    return param_->p;
     	    break;
     	case(NU):
     	    return param_->nu;
     	    break;
     	case(GAMMA):
     	    return param_->gamma;
     	    break;
			case(SIGMA):
					return sigma_;
					break;
     	default:
     	    return -1;
     	    break;
    }
	}

	void SVMWrapper::setTrainingSample(svm_problem* training_sample)
	{
  	training_set_ = training_sample;		
	}
	
	int SVMWrapper::train(struct svm_problem* problem)
	{
	  if (problem != NULL 
				&& param_ != NULL 
				&& (svm_check_parameter(problem,param_) == NULL))
	  {
	  	training_set_ = problem;
	  	
			if (model_ != NULL)
			{
		    svm_destroy_model(model_);
		    model_ = NULL;
			}
			
			if (kernel_type_ == OLIGO)
			{
				if (border_length_ != gauss_table_.size())
				{
					SVMWrapper::calculateGaussTable(border_length_, sigma_, gauss_table_);
    		}
    		training_problem_ = computeKernelMatrix(problem, problem);
    		problem = training_problem_;
			}
			
			model_ = svm_train(problem, param_);
			
			return 1;
		}
	  else
	  {
	  	if (problem == NULL)
	  	{
				cout << "problem is null" << endl;	  		
	  	}
	  	if (param_ == NULL)
	  	{
	  		cout << "param_ == null" << endl;
	  	}
	  	if (svm_check_parameter(problem,param_) != NULL)
	  	{
	  		cout << "check parameter failed: " << endl 
	  		<< svm_check_parameter(problem,param_) << endl;
	  	}
	  	cout << "Training error" << endl;
			return 0;
	  }
	}
	
	void SVMWrapper::saveModel(string model_filename) const throw (Exception::UnableToCreateFile)
	{
		SignedInt status = 0;
		
	  if (model_ != NULL)
	  {
	    status = svm_save_model(model_filename.c_str(), model_);
	  }
	  else
	  {
	  	throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, model_filename);
	  }
	  if (status == -1)
	  {
	  	throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, model_filename);
	  }
	}
	
	void SVMWrapper::loadModel(string model_filename)
	{
		TextFile file;
		TextFile::iterator it;
		vector<String> parts;
		
	  if (model_ != NULL)
	  {
			svm_destroy_model(model_);
			model_ = NULL;
	  }
	  model_ = svm_load_model(model_filename.c_str());
	  setParameter(SVM_TYPE, svm_get_svm_type(model_));
	  file.load(model_filename, true);
	  
		it = file.search("kernel_type");
		if (it != file.end())
		{
			it->split(' ', parts);
			if (parts[1] == "linear")
			{
				setParameter(KERNEL_TYPE, LINEAR);
			}
			else if (parts[1] == "polynomial")
			{
				setParameter(KERNEL_TYPE, POLY);
			}
			else if (parts[1] == "rbf")
			{
				setParameter(KERNEL_TYPE, RBF);
			}
			else if (parts[1] == "sigmoid")
			{
				setParameter(KERNEL_TYPE, SIGMOID);
			}
			else if (parts[1] == "precomputed")
			{
				setParameter(KERNEL_TYPE, OLIGO);
			}
		}
	}
	
	vector<DoubleReal>* SVMWrapper::predict(struct svm_problem* problem)
	{
    DoubleReal          label = 0.0;
    vector<DoubleReal>* results = new vector<DoubleReal>();

		if (model_ == NULL)
			{
				cout << "Model is null" << endl;
			}
		if (problem == NULL)
			{
				cout << "problem is null" << endl;
			}
		if (param_->kernel_type == PRECOMPUTED && training_set_ == NULL)
			{
				cout << "Training set is null and kernel type == PRECOMPUTED" << endl;
			}

    if (model_ == NULL || problem == NULL)
    {
			return results;
    }
		if (kernel_type_ == OLIGO)
		{
			if (training_set_ == NULL)
			{
				return results;
			}
   		problem = computeKernelMatrix(problem, training_set_);
		}
   
    for(int i = 0; i < problem->l; i++)
    {
			label = svm_predict(model_, problem->x[i]);
			results->push_back(label);
		}
		
		if (kernel_type_ == OLIGO)
		{
    	destroyProblem(problem);
		}

    return results;
	}
	
	vector<svm_problem*>* SVMWrapper::createRandomPartitions(svm_problem* problem,
																													 UnsignedInt  number)
	{
		vector<svm_problem*>* problems = new vector<svm_problem*>();
		vector<UnsignedInt> indices;
		UnsignedInt partition_count = 0;
		UnsignedInt actual_partition_size = 0; 
		vector<UnsignedInt>::iterator indices_iterator;
			
		if (number == 1)
		{
			problems->push_back(problem);
			return problems;
		}
		if (number > 1)
		{
			// Creating the particular partition instances
			for(UnsignedInt i = 0; i < number; i++)
			{
				problems->push_back(new svm_problem());
			}
			
			// Creating indices
			for(SignedInt i = 0; i < problem->l; i++)
			{
				indices.push_back(i);
			}
			// Shuffling the indices => random indices
			random_shuffle(indices.begin(), indices.end());
			
			indices_iterator = indices.begin();
			
			for(UnsignedInt partition_index = 0; 
					partition_index < number; 
					partition_index++)
			{
				actual_partition_size = 0;
				// determining the number of elements in this partition
				partition_count = (problem->l / number);
				if (problem->l % number > partition_index)
				{
					partition_count++;
				}
				
				// filling the actual partition with 'partition_count' elements
				while(actual_partition_size < partition_count)
				{
					if (actual_partition_size == 0)
					{
						(*problems)[partition_index]->l = partition_count;
						(*problems)[partition_index]->x = new svm_node*[partition_count];
						(*problems)[partition_index]->y = new DoubleReal[partition_count];						
					}
					(*problems)[partition_index]->x[actual_partition_size] = 
						problem->x[*indices_iterator];
					(*problems)[partition_index]->y[actual_partition_size] = 
						problem->y[*indices_iterator];
					actual_partition_size++;
					indices_iterator++;
				}
			}			
		}
		return problems;		
	}
	
	svm_problem* SVMWrapper::mergePartitions(const vector<svm_problem*>* const problems,
								 															UnsignedInt 											except)
	{
		svm_problem* merged_problem = NULL;
		UnsignedInt count = 0;
		UnsignedInt actual_index = 0;
		
		if (problems->size() == 1 && except == 0)
		{
			return NULL;
		}
		
		if (problems->size() > 0)
		{
			merged_problem = new svm_problem();
			for(UnsignedInt i = 0; i < problems->size(); i++)
			{
				if (i != except)
				{
					count += (*problems)[i]->l;
				}
			}
			merged_problem->l = count;
			merged_problem->x = new svm_node*[count];
			merged_problem->y = new DoubleReal[count];
			for(UnsignedInt i = 0; i < problems->size(); i++)
			{
				if (i != except)
				{
					for(SignedInt j = 0; j < (*problems)[i]->l; j++)
					{
						merged_problem->x[actual_index] = (*problems)[i]->x[j];
						merged_problem->y[actual_index] = (*problems)[i]->y[j];
						actual_index++;
					}
				}
			}
		}
		return merged_problem;
	}
	
	vector<DoubleReal>* SVMWrapper::getLabels(svm_problem* problem)
	{
		UnsignedInt count = 0;
		vector<DoubleReal>* labels = new vector<DoubleReal>();
		
		if (problem == NULL)
		{
			return labels;
		}
		
		count = problem->l;
		for(UnsignedInt i = 0; i < count; i++)
		{
			labels->push_back(problem->y[i]);
		}
		return labels;
	}
	
	map<SVM_parameter_type, DoubleReal>* SVMWrapper::performCrossValidation(svm_problem*   problem,
																 									map<SVM_parameter_type, DoubleReal>&   start_values_map,
																 									map<SVM_parameter_type, DoubleReal>&   step_sizes_map,
																 									map<SVM_parameter_type, DoubleReal>&   end_values_map,
																 									DoubleReal* 												   cv_quality,
																 									UnsignedInt 												   number_of_partitions,
																 									UnsignedInt 												   number_of_runs,
																 									bool																	 additive_step_sizes,
																 									bool				 												   output,
																 									String																 performances_file_name)
	{
		map<SVM_parameter_type, DoubleReal>::iterator start_values_iterator;
		map<SVM_parameter_type, DoubleReal>::iterator step_sizes_iterator;
		map<SVM_parameter_type, DoubleReal>::iterator end_values_iterator;
		map<SVM_parameter_type, DoubleReal>* best_parameters;	
		vector<pair<DoubleReal, UnsignedInt> > combined_parameters;
		combined_parameters.push_back(make_pair(1, 25));
		for(UnsignedInt i = 1; i < gauss_tables_.size(); ++i)
		{
			combined_parameters.push_back(make_pair(1, 25));
		}
		
		DoubleReal precision = 0.0001;
		
		DoubleReal* start_values = new DoubleReal[start_values_map.size()]();
		DoubleReal* actual_values = new DoubleReal[start_values_map.size()]();
		DoubleReal* step_sizes = new DoubleReal[start_values_map.size()]();
		DoubleReal* end_values = new DoubleReal[start_values_map.size()]();
		SVM_parameter_type* actual_types = new SVM_parameter_type[start_values_map.size()]();
		UnsignedInt actual_index = 0;
		bool condition = false;
		bool found = false;
		UnsignedInt counter = 0;
		vector<svm_problem*>* partitions;
		svm_problem** training_data;
		DoubleReal temp_performance = 0;
		vector<DoubleReal>* predicted_labels;
		vector<DoubleReal>* real_labels;
		vector<DoubleReal> performances;
		UnsignedInt max_index = 0;
		DoubleReal max = 0;
		ofstream performances_file;
		ofstream run_performances_file;
		DoubleReal max_performance;
		DoubleReal* best_values = new DoubleReal[start_values_map.size()]();
		vector<DoubleReal>::iterator predicted_it;
		vector<DoubleReal>::iterator real_it;
		
		if (output)
		{
		  performances_file.open(performances_file_name.c_str(), ios_base::out);
		  run_performances_file.open((performances_file_name + "_runs.txt").c_str(), ios_base::out);
		}   
		
		start_values_iterator = start_values_map.begin();
		step_sizes_iterator = step_sizes_map.begin();
		end_values_iterator = end_values_map.begin();

		// Initializing the necessary variables
		while(start_values_iterator != start_values_map.end())
		{
			actual_types[actual_index] = start_values_iterator->first;
			start_values[actual_index] = start_values_iterator->second;
			actual_values[actual_index] = start_values_iterator->second;
			step_sizes_iterator = step_sizes_map.find(start_values_iterator->first);
			if (step_sizes_iterator == step_sizes_map.end())
			{
				
			}
			else
			{
				step_sizes[actual_index] = step_sizes_iterator->second;
			}
			end_values_iterator = end_values_map.find(start_values_iterator->first);
			if (end_values_iterator == end_values_map.end())
			{
				
			}
			else
			{
				end_values[actual_index] = end_values_iterator->second;
			}			
			start_values_iterator++;
			actual_index++;
		}

		// for every 
		for(UnsignedInt i = 0; i < number_of_runs; i++)
		{
			for(UnsignedInt index = 0; index < start_values_map.size(); ++index)
			{
				best_values[index] = 0;
			}
			max_performance = 0;		
			partitions = createRandomPartitions(problem, number_of_partitions);
	
			counter = 0;
			found = true;

			training_data = new svm_problem*[number_of_partitions];
			for(UnsignedInt j = 0; j < number_of_partitions; j++)
			{
				training_data[j] = SVMWrapper::mergePartitions(partitions, j);
			}
			
			while(found)
			{
				// testing whether actual parameters are in the defined range
				condition = true;	
				actual_index = 0;
				while(actual_index < start_values_map.size())
				{			
					if (actual_values[actual_index] > end_values[actual_index])
					{
						condition = false;
					}
					actual_index++;
				}
				// setting the actual parameters
				actual_index = 0;
				while(actual_index < start_values_map.size())
				{
					setParameter(actual_types[actual_index], actual_values[actual_index]);
					
					++actual_index;
					
				}

				// evaluation of parameter performance
				temp_performance = 0;
				for(UnsignedInt j = 0; j < number_of_partitions; j++)
				{
					if (train(training_data[j]))
					{
						predicted_labels = predict((*partitions)[j]);
						real_labels = getLabels((*partitions)[j]);
												
						predicted_it = predicted_labels->begin();
						real_it = real_labels->begin();
						
						if (param_->svm_type == C_SVC || param_->svm_type == NU_SVC)
						{
							temp_performance += Math::BasicStatistics<DoubleReal>::classificationRate(
								predicted_labels->begin(), predicted_labels->end(),
								real_labels->begin(), real_labels->end());
						}
						else if (param_->svm_type == NU_SVR || param_->svm_type == EPSILON_SVR)
						{
							temp_performance += Math::BasicStatistics<DoubleReal>::pearsonCorrelationCoefficient(
								predicted_labels->begin(), predicted_labels->end(),
								real_labels->begin(), real_labels->end());
						}
						
						delete predicted_labels;
						delete real_labels;
						if (param_->kernel_type == PRECOMPUTED)
						{
							destroyProblem(training_problem_);
						}

						if (output && j == number_of_partitions - 1)
						{
							performances_file << temp_performance / (j + 1) << " ";
							for(UnsignedInt k = 0; k < start_values_map.size(); k++)
							{
								switch(actual_types[k])
								{
									case C:
										performances_file << "C: " << actual_values[k];						
										break;
									case NU:
										performances_file << "NU: " << actual_values[k];						
										break;									
									case DEGREE:
										performances_file << "DEGREE: " << actual_values[k];						
										break;																														
									case P:
										performances_file << "P: " << actual_values[k];						
										break;																														
									case GAMMA:										
										performances_file << "GAMMA: " << actual_values[k];						
										break;																														
									case SIGMA:
										performances_file << "SIGMA: " << actual_values[k];																										
										break;																														
									default:
										break;																														
								}
								if (k < (start_values_map.size() - 1))
								{
									performances_file << " ";
								}
								else
								{
									performances_file << endl;
								}
							}
						}
					}
					else
					{
						cout << "Training failed" << endl;
					}
				}

				// storing performance for this parameter combination
				temp_performance = temp_performance / number_of_partitions;
				if (temp_performance > max_performance)
				{
					max_performance = temp_performance;
					for(UnsignedInt index = 0; index < start_values_map.size(); ++index)
					{
						best_values[index] = actual_values[index];
					}		
				}
								
				if (i == 0)
				{
					performances.push_back(temp_performance);
				}
				else
				{
					performances[counter] = performances[counter] + temp_performance;
					counter++;
				}

				// trying to find new parameter combination
				found = false;
				actual_index = 0;
				while(actual_index < start_values_map.size()
							&& !found)
				{
					if (additive_step_sizes)
					{			
						if (actual_values[actual_index] + 
								step_sizes[actual_index] 
								<= end_values[actual_index] + precision)
						{
							found = true;
							actual_values[actual_index] = actual_values[actual_index] + 
																						step_sizes[actual_index];
						}
						else
						{
							actual_values[actual_index] = start_values[actual_index];
						}
					}
					else
					{
						if (actual_values[actual_index] * 
								step_sizes[actual_index] 
								<= end_values[actual_index] + precision)
						{
							found = true;
							actual_values[actual_index] = actual_values[actual_index] * 
																						step_sizes[actual_index];
						}
						else
						{
							actual_values[actual_index] = start_values[actual_index];
						}
					}
					actual_index++;
				}
			}
			
			for(UnsignedInt k = 0; k < number_of_partitions; k++)
			{
				free(training_data[k]->x);
				free(training_data[k]->y);
				delete training_data[k];
			}
			delete training_data;
			if (output)
			{
				cout << "run finished, time elapsed since start: " << clock() 
					<< " mean performance is: " << *(max_element(performances.begin(), performances.end())) / (i + 1) 
					<< endl << "performance of this run is: " << max_performance << " with parameters: ";
				run_performances_file << max_performance << " ";
				for(UnsignedInt k = 0; k < start_values_map.size(); k++)
				{
					switch(actual_types[k])
					{
						case C:
							cout << "C: " << best_values[k];						
							run_performances_file << "C: " << best_values[k];						
							break;
						case NU:
							cout << "NU: " << best_values[k];						
							run_performances_file << "NU: " << best_values[k];						
							break;									
						case DEGREE:
							cout << "DEGREE: " << best_values[k];						
							run_performances_file << "DEGREE: " << best_values[k];						
							break;																														
						case P:
							cout << "P: " << best_values[k];						
							run_performances_file << "P: " << best_values[k];						
							break;																														
						case GAMMA:										
							cout << "GAMMA: " << best_values[k];						
							run_performances_file << "GAMMA: " << best_values[k];						
							break;																														
						case SIGMA:
							cout << "SIGMA: " << best_values[k];						
							run_performances_file << "SIGMA: " << best_values[k];																										
							break;																														
						default:
							break;																														
					}
					if (k < (start_values_map.size() - 1))
					{
						cout << " ";
						run_performances_file << " ";
					}
					else
					{
						cout << endl;
						run_performances_file << endl;
					}
				}

			}
		}
		
		// Determining the index for the maximum performance
		for(UnsignedInt i = 0; i < performances.size(); i++)
		{
			if (performances[i] > max)
			{
				max_index = i;
				max = performances[i];
			}
		}

		// Determining the best parameter combination		
		start_values_iterator = start_values_map.begin();
		actual_index = 0;
		// resetting actual values to start values
		while(start_values_iterator != start_values_map.end())
		{
			actual_values[actual_index] = start_values_iterator->second;
			++actual_index;
			++start_values_iterator;		
		}

		actual_index = 0;
		if (max_index == 0)
		{
			while(actual_index < start_values_map.size())
			{
				actual_values[actual_index] = start_values[actual_index];
				++actual_index;
			}
		}
		else
		{
			performances_file << "Best parameter combination *********************" << endl;
			counter = 1;
			found = true;
			while(found)
			{
				found = false;
				actual_index = 0;
				while(actual_index < start_values_map.size()
							&& !found && counter <= max_index)
				{			
					if (additive_step_sizes)
					{			
						if (actual_values[actual_index] + 
								step_sizes[actual_index] 
								<= end_values[actual_index] + precision)
						{
							found = true;
							actual_values[actual_index] = actual_values[actual_index] + 
																						step_sizes[actual_index];
						}
						else
						{
							actual_values[actual_index] = start_values[actual_index];
						}
					}
					else
					{
						if (actual_values[actual_index] * 
								step_sizes[actual_index] 
								<= end_values[actual_index] + precision)
						{
							found = true;
							actual_values[actual_index] = actual_values[actual_index] * 
																						step_sizes[actual_index];
						}
						else
						{
							actual_values[actual_index] = start_values[actual_index];
						}
					}
					actual_index++;
				}
				performances_file	<< performances[counter]  / number_of_runs << ": ";
				for(UnsignedInt k = 0; k < start_values_map.size(); k++)
				{
					switch(actual_types[k])
					{
						case C:
							performances_file << "C: " << actual_values[k];						
							break;
						case NU:
							performances_file << "NU: " << actual_values[k];						
							break;									
						case DEGREE:
							performances_file << "DEGREE: " << actual_values[k];						
							break;																														
						case P:
							performances_file << "P: " << actual_values[k];						
							break;																														
						case GAMMA:										
							performances_file << "GAMMA: " << actual_values[k];						
							break;																														
						case SIGMA:
							performances_file << "SIGMA: " << actual_values[k];																										
							break;																														
						default:
							break;
					}																														
					if (k < (start_values_map.size() - 1))
					{
						performances_file << ", ";
					}
					else
					{
						performances_file << endl;
					}
				}
				// best parameter combination found
				if (counter == max_index)
				{
					found = false;
				}
				counter++;
			}
			performances_file << "Best parameter combination ended****************" << endl;
		}
		
		performances_file.close();
		
		best_parameters = new map<SVM_parameter_type, DoubleReal>();
		
		for(actual_index = 0; actual_index < start_values_map.size(); actual_index++)
		{
			best_parameters->insert(make_pair(actual_types[actual_index], actual_values[actual_index]));
		}
		*cv_quality = performances[max_index] / number_of_runs;

		actual_index = 0;
		while(actual_index < start_values_map.size())
		{			
			setParameter(actual_types[actual_index], actual_values[actual_index]);
			actual_index++;
		}
				
		delete start_values;
		delete actual_values;
		delete end_values;
		delete step_sizes;
		delete best_values;
		
		return best_parameters;									
	}

	vector<DoubleReal>* SVMWrapper::predict(const std::vector<svm_node*>& vectors)
	{
    DoubleReal          label = 0.0;
    vector<DoubleReal>* results = new vector<DoubleReal>();

    if (model_ == NULL)
    {
			return results;
    }
   
    for(UnsignedInt i = 0; i < vectors.size(); i++)
    {
			label = svm_predict(model_, vectors[i]);
			results->push_back(label);
		}
		
    return results;
	}
	
	double SVMWrapper::getSVRProbability()
	{
		if (model_ != NULL)
		{
			return svm_get_svr_probability(model_);
		}
		else
		{
			return 0;
		}
	}																			 						

	void SVMWrapper::initParameters()
	{
		model_ = NULL;
	  
	  param_->svm_type = NU_SVR;
	  param_->kernel_type = PRECOMPUTED;
	  param_->degree = 1;                            // for poly
	  param_->gamma = 1.0;	                         // for poly/rbf/sigmoid 
	  param_->coef0 = 0;	                           // for poly/sigmoid 
	  param_->cache_size = 300;                      // in MB 
	  param_->eps = 0.001;	                         // stopping criterium 
	  param_->C = 1;	                               // for C_SVC, EPSILON_SVR, and NU_SVR 
	  param_->nu = 0.5;	                             // for NU_SVC, ONE_CLASS, and NU_SVR 
	  param_->p = 0.1;	                             // for EPSILON_SVR 
	  param_->shrinking = 0;	                       // use the shrinking heuristics 
		param_->probability = 0;

		param_->nr_weight = 0;
	}

	DoubleReal SVMWrapper::kernelOligo(const svm_node* 						x, 
																		 const svm_node*						y,
																		 const vector<DoubleReal>& 	gauss_table,
																		 DoubleReal 								sigma_square,
												  					 UnsignedInt 								max_distance)
  {
    double kernel = 0;
    int    i1     = 0;
    int    i2     = 0;
    int    c1     = 0;

    while(x[i1].index != -1
	&& y[i2].index != -1)
    {
      if (x[i1].index == y[i2].index)
      {
  	if (((UnsignedInt) abs(x[i1].value - y[i2].value)) <= max_distance)
    	{
          if (sigma_square == 0)
          {
            try
            {
              kernel += gauss_table.at(abs((int)(x[i1].value
                                                - y[i2].value)));
            }
            catch(...)
            {
              cout << "Tried to access " << x[i1].value << " - " << y[i2].value << endl;            
            }
          }
          else
          {
//						cout << "kernelOligo" << endl;
            kernel += exp(-1 * (x[i1].value - y[i2].value)* (x[i1].value - y[i2].value) / (4 * sigma_square));

//            cout << "adding " << exp(-1 * (x[i1].value - y[i2].value)* (x[i1].value - y[i2].value) / (4 * sigma_square))
//              << " for index " << x[i1].index << " and positions " << x[i1].value << ", " << y[i2].value << endl;

          }
          if (x[i1].index == x[i1 + 1].index)
          {
            i1++;
            c1++;
          }
          else if (y[i2].index == y[i2 + 1].index)
          {
            i2++;
            i1 -= c1;
            c1 = 0;
          }
          else
          {
            i1++;
            i2++;
          }
   	}
   	else
   	{
          if (x[i1].value < y[i2].value)
          {
            if (x[i1].index == x[i1 + 1].index)
            {
              i1++;
            }
            else if (y[i2].index == y[i2 + 1].index)
            {
              i2++;
              i1 -= c1;
              c1 = 0;
            }
            else
            {
              i1++;
              i2++;
            }
          }
          else
          {
            i2++;
            i1 -= c1;
            c1 = 0;
          }
        }
      }
      else
      {
        if (x[i1].index < y[i2].index)
        {
          i1++;
        }
        else
        {
          i2++;
        }
        c1 = 0;
      }
    }
    return kernel;
  }		

	svm_problem* SVMWrapper::computeKernelMatrix(svm_problem* problem1, svm_problem* problem2)
	{
		DoubleReal temp = 0;
		svm_problem* kernel_matrix;
		vector<DoubleReal> sigma_squares;
				
		if (problem1 == NULL || problem2 == NULL)
		{
			return NULL;
		}	
		UnsignedInt number_of_sequences = 0;

		number_of_sequences = problem1->l;		
		kernel_matrix = new svm_problem;
		kernel_matrix->l = number_of_sequences;
		kernel_matrix->x = new svm_node*[number_of_sequences];
		kernel_matrix->y = new DoubleReal[number_of_sequences];
		
		for(UnsignedInt i = 0; i < number_of_sequences; i++)
		{
			kernel_matrix->x[i] = new svm_node[problem2->l + 2];
			kernel_matrix->x[i][0].index = 0;
			kernel_matrix->x[i][0].value = i + 1;
			kernel_matrix->y[i] = problem1->y[i];
			kernel_matrix->x[i][problem2->l + 1].index = -1;
		}

		if (problem1 == problem2)
		{
			for(UnsignedInt i = 0; i < number_of_sequences; i++)
			{			
				for(UnsignedInt j = i; j < number_of_sequences; j++)
				{
					temp = SVMWrapper::kernelOligo(problem1->x[i], problem2->x[j], gauss_table_);
					kernel_matrix->x[i][j + 1].index = j + 1;
					kernel_matrix->x[i][j + 1].value = temp;
					kernel_matrix->x[j][i + 1].index = i + 1;
					kernel_matrix->x[j][i + 1].value = temp;				
				}
			}
		}
		else
		{
			for(UnsignedInt i = 0; i < number_of_sequences; i++)
			{			
				for(UnsignedInt j = 0; j < (UnsignedInt) problem2->l; j++)
				{
					temp = SVMWrapper::kernelOligo(problem1->x[i], problem2->x[j], gauss_table_);

					kernel_matrix->x[i][j + 1].index = j + 1;
					kernel_matrix->x[i][j + 1].value = temp;
				}
			}			
		}
		return kernel_matrix;
	}
	
	void SVMWrapper::destroyProblem(svm_problem* problem)
	{
		for(SignedInt i = 0; i < problem->l; i++)
		{
			free(problem->x[i]);
		}
		free(problem->y);
		free(problem->x);
		free(problem);
	}
		
	void SVMWrapper::getSignificanceBorders(svm_problem* data, 
																					pair<DoubleReal, DoubleReal>& sigmas,
																					DoubleReal confidence,
																					UnsignedInt number_of_runs,
																					UnsignedInt number_of_partitions,
																					DoubleReal step_size,
																					UnsignedInt max_iterations)
	{
		vector<pair<DoubleReal, DoubleReal> > points;
		vector<DoubleReal> 										differences;
		vector<svm_problem*>* 								partitions;
		svm_problem*													training_data;
		vector<DoubleReal>*										predicted_labels;
		vector<DoubleReal>*										real_labels;
		UnsignedInt														counter = 0;
		UnsignedInt														target = 0;
		ofstream															file("points.txt");
		DoubleReal 														mean;
		DoubleReal														sigma1 = 0;
		DoubleReal														sigma2 = 0;
			
		
		// creation of points (measured rt, predicted rt)
		for(UnsignedInt i = 0; i < number_of_runs; ++i)
		{
			partitions = createRandomPartitions(data, number_of_partitions);
			
			for (UnsignedInt j = 0; j < number_of_partitions; ++j)
			{
				training_data = SVMWrapper::mergePartitions(partitions, j);
				if (train(training_data))
				{
					predicted_labels = predict((*partitions)[j]);
					real_labels = getLabels((*partitions)[j]);
					vector<DoubleReal>::iterator pred_it = predicted_labels->begin();
					vector<DoubleReal>::iterator real_it = real_labels->begin();						
					while(pred_it != predicted_labels->end()
								&& real_it != real_labels->end())
					{
						points.push_back(make_pair(*real_it, *pred_it));
						differences.push_back(abs(*real_it - *pred_it));
						file << *real_it << " " << *pred_it << endl;
						++pred_it;
						++real_it;
					}
					delete predicted_labels;
					delete real_labels;
				}
			}
		}
		file << flush;
								
		// trying to find the two line parameters
		target = (UnsignedInt) round(confidence * points.size());
		
		mean = accumulate(differences.begin(), differences.end(), 0.0) / differences.size();
		sigma1 = mean;
		sigma2 = mean;
		while(target > getNumberOfEnclosedPoints(sigma1, sigma2, points) && counter < max_iterations)
		{
			
			cout << "sigma1: " << sigma1 << ", sigma2: " << sigma2 << " shape contains " 
						<< ((getNumberOfEnclosedPoints(sigma1, sigma2, points) / ((DoubleReal)points.size())) * 100)
						<< " % of points" << endl;

			sigma1 += (step_size / 2);
			sigma2 += step_size;
			++counter;			
		}
		sigmas.first = sigma1;
		sigmas.second = sigma2;

		cout << "sigma1: " << sigma1 << ", sigma2: " << sigma2 << " shape contains " 
			<< ((getNumberOfEnclosedPoints(sigma1, sigma2, points) / ((DoubleReal)points.size())) * 100)
			<< " % of points" << endl;			
	}
	
	UnsignedInt SVMWrapper::getNumberOfEnclosedPoints(DoubleReal sigma1, 
																										DoubleReal sigma2, 
																										const vector<pair<DoubleReal, DoubleReal> >& points)
	{
		UnsignedInt counter = 0;
		DoubleReal 	sigma		= 0;
		
		for(vector<pair<DoubleReal, DoubleReal> >::const_iterator it = points.begin();
				it != points.end();
				++it)
		{
			sigma = sigma1 + it->first * (sigma2 - sigma1);
			
			if (it->first + sigma >= it->second
					&& it->first - sigma <= it->second)
			{
				++counter;
			}
		}
		return counter; 
	}
	
	DoubleReal SVMWrapper::getPValue(DoubleReal 										sigma1, 
																	 DoubleReal 										sigma2,
																	 pair<DoubleReal, DoubleReal>		point)
	{
		DoubleReal center = point.first;
		DoubleReal distance = point.second - center;
		DoubleReal sd_units = 0;
		DoubleReal result = 0;
		
		// getting the absolute value
		distance = (distance < 0) ? (-1 * distance) : distance;
			
		sd_units = (2 * distance) / (sigma1 + point.first * (sigma2 - sigma1));
		
		result = (gsl_cdf_gaussian_P(sd_units, 1) - 0.5) * 2;
		if (result > 1)
		{
			result = 1;
		}
		return result;
	}
	
	void SVMWrapper::getDecisionValues(svm_problem* data, vector<DoubleReal>& decision_values)
	{
		DoubleReal temp_value;
		
		decision_values.clear();
		if (model_ != NULL)
		{
			for(SignedInt i = 0; i < data->l; ++i)
			{
				temp_value = 0;						
				svm_predict_values(model_, data->x[i], &temp_value);
				decision_values.push_back(temp_value);
			}
		}
	}																		  			
	
	void SVMWrapper::scaleData(svm_problem* data, SignedInt max_scale_value)
	{
		vector<DoubleReal> max_values;
		vector<DoubleReal> min_values;
		vector<DoubleReal> sums;
		SignedInt max_index = 0;
		SignedInt j = 0;
		
		for(SignedInt i = 0; i < data->l; ++i)
		{
			j = 0;
			while(data->x[i][j].index != -1)
			{
				if (data->x[i][j].index > max_index)
				{
					max_index = data->x[i][j].index;
				}
				++j;
			}
		}
		
		max_values.resize(max_index, 0);
		min_values.resize(max_index, 0);
		sums.resize(max_index, 0);
		
		for(SignedInt i = 0; i < data->l; ++i)
		{
			j = 0;
			while(data->x[i][j].index != -1)
			{
				if (data->x[i][j].value > max_values.at(data->x[i][j].index - 1))
				{
					max_values.at(data->x[i][j].index - 1) = data->x[i][j].value;
				}
				sums.at(data->x[i][j].index - 1) = sums.at(data->x[i][j].index - 1) + data->x[i][j].value;
				if (data->x[i][j].value < min_values.at(data->x[i][j].index - 1))
				{
					min_values.at(data->x[i][j].index - 1) = data->x[i][j].value;
				}

				++j;
			}
		}
		for(SignedInt i = 0; i < data->l; ++i)
		{
			j = 0;
			while(data->x[i][j].index != -1)				
			{
				if (max_scale_value == -1)
				{
					data->x[i][j].value = 2 * (data->x[i][j].value - min_values.at(data->x[i][j].index - 1))
																/ (max_values.at(data->x[i][j].index - 1) - min_values.at(data->x[i][j].index - 1)) - 1;																			
				}
				else
				{
					data->x[i][j].value = max_scale_value * (data->x[i][j].value - min_values.at(data->x[i][j].index - 1))
																	/ (max_values.at(data->x[i][j].index - 1) - min_values.at(data->x[i][j].index - 1));																			
				}
				++j;
			}
		}
	}																		  																										

	void SVMWrapper::calculateGaussTable(UnsignedInt border_length, 
																			 DoubleReal sigma, 
																			 vector<DoubleReal>&	gauss_table)
	{
		if (border_length != gauss_table.size())
		{
			gauss_table.resize(border_length, 0);
		}		
		gauss_table[0] = 1;
	 	for(UnsignedInt i = 1; i < border_length; ++i)
	 	{
	  	gauss_table[i] = exp((-1 / 4.0 /
					 						     (sigma * sigma)) *
	  								 	     (i * i));
	  }
	}
	

} // namespace OpenMS
