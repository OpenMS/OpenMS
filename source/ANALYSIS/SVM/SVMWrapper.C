// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

namespace OpenMS 
{
		
	SVMWrapper::SVMWrapper()
	{
	
	  model_ = NULL;
	  param_ = (struct svm_parameter*) malloc(sizeof(struct svm_parameter));
	  initParameters();
	}
	
	SVMWrapper::~SVMWrapper()
	{
	  if (param_ != NULL)
	  {
			free(param_);                           
	  }
	  if (model_ != NULL)
	  {
			svm_destroy_model(model_);
	  } 
	}
	
	void SVMWrapper::setParameter(SVM_parameter_type type, int value)
	{
	
	  switch(type)
	  {
	  	case(SVM_TYPE):
				if (value == NU_SVR || value == EPSILON_SVR)	
				{
				  param_->svm_type = value;
				}
				break;
			case(KERNEL_TYPE):	
				param_->kernel_type = value;
				break;
			case(DEGREE):
				param_->degree = value;
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
	    	    return param_->kernel_type;
	    	    break;
	    	case(SVM_TYPE):
	    	    return param_->svm_type;
	    	    break;
	    	case(DEGREE):
	    	    return (int) param_->degree;
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
				case(C):
			    param_->C = value;
		   		break;
				case(P):
			    param_->p = value;
		   		break;
				case(NU):
			    param_->nu = value;
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
     	default:
     	    return -1;
     	    break;
    }
	}
	
	int SVMWrapper::train(struct svm_problem* problem)
	{
	  if (problem != NULL 
				&& param_ != NULL 
				&& (svm_check_parameter(problem,param_) == NULL))
	  {
			if (model_ != NULL)
			{
		    svm_destroy_model(model_);
		    model_ = NULL;
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
	
	void SVMWrapper::saveModel(string model_filename)
	{
	  if (model_ != NULL)
	  {
	    svm_save_model(model_filename.c_str(), model_);
	  }
	}
	
	void SVMWrapper::loadModel(string model_filename)
	{
	  if (model_ != NULL)
	  {
			svm_destroy_model(model_);
			model_ = NULL;
	  }
	  model_ = svm_load_model(model_filename.c_str());
	}
	
	vector<DoubleReal>* SVMWrapper::predict(struct svm_problem* problem)
	{
    DoubleReal          label = 0.0;
    vector<DoubleReal>* results = new vector<DoubleReal>();

    if (model_ == NULL || problem == NULL)
    {
			return results;
    }
   
    for(int i = 0; i < problem->l; i++)
    {
			label = svm_predict(model_, problem->x[i]);
			results->push_back(label);
		}
		
    return results;
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
																 									bool				 												   output)
	{
		map<SVM_parameter_type, DoubleReal>::iterator start_values_iterator;
		map<SVM_parameter_type, DoubleReal>::iterator step_sizes_iterator;
		map<SVM_parameter_type, DoubleReal>::iterator end_values_iterator;
		map<SVM_parameter_type, DoubleReal>* best_parameters;	
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
		ofstream performances_file("performances.txt");
		
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
					actual_index++;
				}

				// evaluation of parameter performance
				temp_performance = 0;
				for(UnsignedInt j = 0; j < number_of_partitions; j++)
				{
					if (train(training_data[j]))
					{
						predicted_labels = predict((*partitions)[j]);
						real_labels = getLabels((*partitions)[j]);
						vector<DoubleReal>::iterator predicted_it = predicted_labels->begin();
						vector<DoubleReal>::iterator real_it = real_labels->begin();

						temp_performance += Math::BasicStatistics<DoubleReal>::pearsonCorrelationCoefficient(
							predicted_labels->begin(), predicted_labels->end(),
							real_labels->begin(), real_labels->end());
						
						delete predicted_labels;
						delete real_labels;												
					}
				}

				// storing performance for this parameter combination
				temp_performance = temp_performance / number_of_partitions;				
				if (i == 0)
				{
					performances.push_back(temp_performance);
				}
				else
				{
					performances[counter] = performances[counter] + temp_performance;
					counter++;
				}

				// trying to set new parameter combination
				found = false;
				actual_index = 0;
				while(actual_index < start_values_map.size()
							&& !found)
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
				cout << "run finished, time elapsed since start: " << clock() << endl;
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
			actual_index++;
			start_values_iterator++;		
		}

		actual_index = 0;
		if (max_index == 0)
		{
			while(actual_index < start_values_map.size())
			{
				actual_values[actual_index] = start_values[actual_index];
				actual_index++;
			}
		}
		else
		{
			counter = 1;
			found = true;
			while(found)
			{
				found = false;
				actual_index = 0;
				while(actual_index < start_values_map.size()
							&& !found && counter <= max_index)
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
					actual_index++;
				}
				performances_file << performances[counter]  / number_of_runs << ": ";
				for(UnsignedInt k = 0; k < start_values_map.size(); k++)
				{
					if (actual_types[k] == C)
					{
						performances_file << "C: " << actual_values[k];						
					}
					else if (actual_types[k] == NU)
					{
						performances_file << "NU: " << actual_values[k];						
					}
					else if (actual_types[k] == DEGREE)
					{
						performances_file << "DEGREE: " << actual_values[k];						
					}
					else if (actual_types[k] == P)
					{
						performances_file << "P: " << actual_values[k];						
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
		
		return best_parameters;									
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
	  param_->kernel_type = POLY;
	  param_->degree = 1;                            // for poly
	  param_->gamma = 1.0;	                         // for poly/rbf/sigmoid 
	  param_->coef0 = 0;	                           // for poly/sigmoid 
	  param_->cache_size = 300;                      // in MB 
	  param_->eps = 0.001;	                         // stopping criterium 
	  param_->C = 1;	                               // for C_SVC, EPSILON_SVR, and NU_SVR 
	  param_->nr_weight = 0;	                       // for C_SVC 
	  param_->nu = 0.5;	                             // for NU_SVC, ONE_CLASS, and NU_SVR 
	  param_->p = 0.1;	                             // for EPSILON_SVR 
	  param_->shrinking = 0;	                       // use the shrinking heuristics 
		param_->probability = 0;
	}

} // namespace OpenMS
