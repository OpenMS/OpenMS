// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Sandro Andreotti $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_SVM_SVMWRAPPER_H
#define OPENMS_ANALYSIS_SVM_SVMWRAPPER_H

#include <svm.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>

namespace OpenMS 
{
	
	/// Data structure used in SVMWrapper
	struct SVMData
	{
		std::vector< std::vector< std::pair<Int, DoubleReal> > > sequences;
		std::vector< DoubleReal > labels;
		
		SVMData()
			: sequences(std::vector< std::vector< std::pair<Int, DoubleReal> > >()),
				labels(std::vector< DoubleReal >())
		{
		};
		
		SVMData(std::vector< std::vector< std::pair<Int, DoubleReal> > >& seqs, std::vector< DoubleReal >& lbls)
			: sequences(seqs),
				labels(lbls)
		{
		};
		
		bool operator == (const SVMData& rhs) const
		{
			return sequences == rhs.sequences
				&& labels == rhs.labels;
		}
		
		bool store(const String& filename) const
		{
			std::ofstream output_file(filename.c_str());

			// checking if file is writable
			if (!File::writable(filename) || sequences.size() != labels.size())
			{
				return false;
			}
				
			// writing feature vectors		
			for (Size i = 0; i < sequences.size(); i++)
			{
				output_file << labels[i] << " ";
				for (Size j = 0; j < sequences[i].size(); ++j)
				{
					output_file << sequences[i][j].second << ":" << sequences[i][j].first << " " ;					
				}
				output_file << std::endl;
			}
			output_file.flush();
			output_file.close();
			std::cout.flush();
			return true;
		}
		
		bool load(const String& filename)
		{
			Size counter = 0;
			std::vector<String> parts;
			std::vector<String> temp_parts;
			
			if (!File::exists(filename))
			{
				return false;
			}
			if (!File::readable(filename))
			{
				return false;
			}
	    if (File::empty(filename))
	    {
				return false;
	    }		
	
			TextFile text_file(filename.c_str(), true);
	    TextFile::iterator it;
	
			it = text_file.begin();
			
			sequences.resize(text_file.size(), std::vector< std::pair<Int, DoubleReal> >());
			labels.resize(text_file.size(), 0.);
			while(counter < text_file.size()&& it != text_file.end())
			{
				it->split(' ', parts);
				labels[counter] = parts[0].trim().toFloat();		
				sequences[counter].resize(parts.size(), std::pair<Int, DoubleReal>());			
				for (Size j = 1; j < parts.size(); ++j)
				{
					parts[j].split(':', temp_parts);
					if (temp_parts.size() < 2)
					{
						return false;
					}
					sequences[counter][j - 1].second = temp_parts[0].trim().toFloat();
					sequences[counter][j - 1].first = temp_parts[1].trim().toInt();
				}
				++counter;
				++it;
			}
			return true;			
		}
							
	};

  /**
    @brief Serves as a wrapper for the libsvm
    
    This class can be used for svm predictions. You can either perform classification or regression and
		choose certain kernel fuctions and additional parameters. Furthermore the models can be saved and 
		loaded and we support also a new kernel function that was specially designed for learning with
		small sequences of different lengths.
  */
  class OPENMS_DLLAPI SVMWrapper
    : public ProgressLogger
	{
	 public:
	
		  /**
		    @brief 	Parameters for the svm to be set from outside
		
				This type is used to specify the kind of parameter that
				is to be set or retrieved by the set/getParameter methods.
			*/
			enum SVM_parameter_type
			{
				SVM_TYPE,				///< the svm type cab be NU_SVR or EPSILON_SVR
				KERNEL_TYPE,		///< the kernel type
				DEGREE,					///< the degree for the polynomial- kernel
				C,							///< the C parameter of the svm
				NU,							///< the nu parameter for nu-SVR
				P,							///< the epsilon parameter for epsilon-SVR
				GAMMA,					///< the gamma parameter of the POLY, RBF and SIGMOID kernel
				PROBABILITY,		///<
				SIGMA,					///<
				BORDER_LENGTH		///<
			};
			
			/// Kernel type
			enum SVM_kernel_type
			{
				OLIGO = 19,
				OLIGO_COMBINED
			};
	
			/// standard constructor
	    SVMWrapper();

			/// destructor	
	    virtual ~SVMWrapper();
	
		  /**
		    @brief You can set the parameters of the svm: 
	     
	         KERNEL_TYPE: can be LINEAR              		 for the linear kernel
	                             RBF                 		 for the rbf kernel
	                             POLY                    for the polynomial kernel
	                             SIGMOID           			 for the sigmoid kernel 
	         DEGREE:      the degree for the polynomial- kernel and the
	                      locality- improved kernel
	     
	         C:            the C parameter of the svm	     
			*/	     
	    void setParameter(SVM_parameter_type type, Int value);
	
		  /**
		    @brief sets the double parameters of the svm
		    
			*/
	    void setParameter(SVM_parameter_type type, DoubleReal value);
		
		  /**
		    @brief	trains the svm 

	      The svm is trained with the data stored in the 'svm_problem' structure.
			*/
	    Int train(struct svm_problem* problem);
	
		  /**
		    @brief	trains the svm 

	      The svm is trained with the data stored in the 'SVMData' structure.
			*/
			Int train(SVMData& problem);

		  /**
		    @brief	saves the svm model 

	      The model of the trained svm is saved into 'modelFilename'. Throws an exception if 
	      the model cannot be saved.

				@exception Exception::UnableToCreateFile
			*/
	    void saveModel(std::string modelFilename) const;
	
		  /**
		    @brief loads the model

	      The svm- model is loaded. After this, the svm is ready for 
	      prediction.
		  */
	    void loadModel(std::string modelFilename);
	
		  /**
		    @brief predicts the labels using the trained model
		    
	     	 The prediction process is started and the results are stored in 'predicted_labels'.

		  */
	    void predict(struct svm_problem* problem, std::vector<DoubleReal>& predicted_labels);

		  /**
		    @brief predicts the labels using the trained model
		    
	     	 The prediction process is started and the results are stored in 'predicted_labels'.

		  */
			void predict(const SVMData& problem, std::vector<DoubleReal>& results);

		  /**
		    @brief You can get the actual int- parameters of the svm

	         KERNEL_TYPE: can be LINEAR              		 for the linear kernel
	                             RBF                 		 for the rbf kernel
	                             POLY                     for the polynomial kernel
	                             SIGMOID           			 for the sigmoid kernel 
	     
	         DEGREE:       the degree for the polynomial- kernel and the
	                       locality- improved kernel
	     
	         SVM_TYPE:     the SVm type of the svm: can be NU_SVR or EPSILON_SVR
		  */	     
	    Int getIntParameter(SVM_parameter_type type);
	
		  /**
		    @brief You can get the actual double- parameters of the svm
		    
		    C:            the C parameter of the svm
		    P:			      the P parameter of the svm (sets the epsilon in
	     											  epsilon-svr)
	   		NU:           the nu parameter in nu-SVR
	   		GAMMA:				for POLY, RBF and SIGMOID		    
		  */
	    DoubleReal getDoubleParameter(SVM_parameter_type type); 

		  /**
		    @brief You can create 'number' equally sized random partitions
		    
		    This function creates 'number' equally sized random partitions and stores them in 'partitions'. 
		    
		  */
			static void createRandomPartitions(svm_problem* problem, Size number, std::vector<svm_problem*>& partitions);
	
		  /**
		    @brief You can create 'number' equally sized random partitions
		    
		    This function creates 'number' equally sized random partitions and stores them in 'partitions'. 
		    
		  */
			static void createRandomPartitions(const SVMData&				  problem,
																				Size  			 					  number,
																				std::vector<SVMData>&		problems);
		  /**
		    @brief You can merge partitions excuding the partition with index 'except' 
		    
		  */
	    static svm_problem* mergePartitions(const std::vector<svm_problem*>& problems, Size except);
																	 				
		  /**
		    @brief You can merge partitions excuding the partition with index 'except' 
		    
		  */
			static void mergePartitions(const std::vector<SVMData>& problems,
			 													  Size 											  except,
				 													SVMData&										merged_problem);

		  /**
		    @brief predicts the labels using the trained model
		    
	     	 The prediction process is started and the results are stored in 'predicted_rts'.

		  */
	    void predict(const std::vector<svm_node*>& vectors, std::vector<DoubleReal>& predicted_rts);
	
		  /**
		    @brief Stores the stored labels of the encoded SVM data at 'labels' 
		    
		  */
			static void getLabels(svm_problem* problem, std::vector<DoubleReal>& labels);
																	 				
		  /**
		    @brief Performs a CV for the data given by 'problem'
		    
		  */
      DoubleReal performCrossValidation(svm_problem*   															      problem_ul,
			            									    const SVMData&		                                problem_l,
									                    	const bool                                        is_labeled,
																 				const	std::map<SVM_parameter_type, DoubleReal>&   start_values_map,
																 				const	std::map<SVM_parameter_type, DoubleReal>&   step_sizes_map,
																 				const	std::map<SVM_parameter_type, DoubleReal>&   end_values_map,
																 				Size     												   				 		    number_of_partitions,
																 				Size 												   				 			      number_of_runs,
																 				std::map<SVM_parameter_type, DoubleReal>&  	      best_parameters,
												 								bool																	 			 			additive_step_sizes = true,
												 								bool				 												   			 			output = false,
												 								String																 			 			performances_file_name = "performances.txt",
												 								bool																				 			mcc_as_performance_measure = false);

      
      /**
		    @brief Returns the probability parameter sigma of the fitted laplace model.		      
		    
		    The libsvm is used to fit a laplace model to the prediction values by performing
		    an internal cv using the training set if setParameter(PROBABILITY, 1) was invoked
		    before using train. Look for your libsvm documentation for more details.
		    The model parameter sigma is returned by this method.	If no model was fitted during 
		    training zero is returned. 
		  */
			DoubleReal getSVRProbability();																			 					

      /**
				@brief returns the value of the oligo kernel for sequences 'x' and 'y'

	        This function computes the kernel value of the oligo kernel,
					which was introduced by Meinicke et al. in 2004. 'x' and
					'y' are encoded by encodeOligo and 'gauss_table' has to be 
					constructed by calculateGaussTable. 
					
					'max_distance' can be used to speed up the computation 
					even further by restricting the maximum distance between a k_mer at
					position i in sequence 'x' and a k_mer at position j 
					in sequence 'y'. If i - j > 'max_distance' the value is not
					added to the kernel value. This approximation is switched
					off by default (max_distance < 0).
      */
      static DoubleReal kernelOligo(const std::vector< std::pair<int, double> >&    x, 
																		const std::vector< std::pair<int, double> >&    y,
																		const std::vector<double>& 	                    gauss_table,
																		int 			                    									max_distance = -1);
																
		  /**
		    @brief calculates the oligo kernel value for the encoded sequences 'x' and 'y'
		    
		    This kernel function calculates the oligo kernel value [Meinicke 04] for
		    the sequences 'x' and 'y' that had been encoded by the encodeOligoBorder... function
		    of the LibSVMEncoder class. 	      		    
		  */
			static DoubleReal kernelOligo(const svm_node*	x, const svm_node* y, const std::vector<DoubleReal>&	gauss_table, DoubleReal sigma_square = 0,	Size	max_distance = 50);

		  /**
		    @brief calculates the significance borders of the error model and stores them in 'sigmas'	      
		    
		  */
			void getSignificanceBorders(svm_problem* data, std::pair<DoubleReal, DoubleReal>& borders, DoubleReal confidence = 0.95, Size number_of_runs = 5, Size number_of_partitions = 5, DoubleReal step_size = 0.01, Size max_iterations = 1000000);

		  /**
		    @brief calculates the significance borders of the error model and stores them in 'sigmas'	      
		    
		  */
			void getSignificanceBorders(const SVMData& data, 
																	std::pair<DoubleReal, DoubleReal>& sigmas,
																	DoubleReal confidence = 0.95,
																	Size number_of_runs = 5,
																	Size number_of_partitions = 5,
																	DoubleReal step_size = 0.01,
																	Size max_iterations = 1000000);

		  /**
		    @brief calculates a p-value for a given data point using the model parameters
		    
		    Uses the model parameters to calculate the p-value for 'point' which has the data
		    entries: measured, predicted retention time.	      
		    
		  */
			DoubleReal getPValue(DoubleReal sigma1, DoubleReal sigma2, std::pair<DoubleReal, DoubleReal> point);

		  /**
		    @brief stores the prediction values for the encoded data in 'decision_values'
		    
		    This function can be used to get the prediction values of the data if a model 
		    is already trained by the train() method. For regression the result is the same
		    as for the method predict. For classification this function returns the distance from
		    the separating hyperplane. For multiclass classification the decision_values vector
		    will be empty. 
		    
		  */
			void getDecisionValues(svm_problem* data, std::vector<DoubleReal>& decision_values);

		  /**
		    @brief Scales the data such that every coloumn is scaled to [-1, 1].
		    
		    Scales the x[][].value values of the svm_problem* structure. If the second 
		    parameter is omitted, the data is scaled to [-1, 1]. Otherwise the data is scaled to [0, max_scale_value]
		  */
			void scaleData(svm_problem* data, Int max_scale_value = -1);

			static void calculateGaussTable(Size border_length, DoubleReal sigma, std::vector<DoubleReal>&	gauss_table);

		  /**
		    @brief computes the kernel matrix using the actual svm parameters	and the given data
		    
		    This function can be used to compute a kernel matrix. 'problem1' and 'problem2'
		    are used together wit the oligo kernel function (could be extended if you 
		    want to use your own kernel functions).      
		    
		  */
			svm_problem* computeKernelMatrix(svm_problem* problem1, svm_problem* problem2);

		  /**
		    @brief computes the kernel matrix using the actual svm parameters	and the given data
		    
		    This function can be used to compute a kernel matrix. 'problem1' and 'problem2'
		    are used together wit the oligo kernel function (could be extended if you 
		    want to use your own kernel functions).      
		    
		  */
			svm_problem* computeKernelMatrix(const SVMData& problem1, const SVMData& problem2);

		  /**
		    @brief This is used for being able to perform predictions with non libsvm standard kernels	    		        
		    
		  */
  		void setTrainingSample(svm_problem* training_sample);
  		
		  /**
		    @brief This is used for being able to perform predictions with non libsvm standard kernels	    		        
		    
		  */
			void setTrainingSample(SVMData& training_sample);
				
		  /**
		    @brief This function fills probabilities with the probability estimates for the first class.
		    		    		        
		    The libSVM function svm_predict_probability is called to get probability estimates
		    for the positive class. Since this is only used for binary classification it is sufficient
		    for every test example to report the probability of the test example belonging to the positive
		    class. Probability estimates have to be turned on during training (svm.setParameter(PROBABILITY, 1)),
		    otherwise this method will fill the 'probabilities' vector with -1s.
		  */
 			void getSVCProbabilities(struct svm_problem* problem, std::vector<DoubleReal>& probabilities, std::vector<DoubleReal>& prediction_labels);

		  /**
		    @brief Sets weights for the classes in C_SVC (see libsvm documentation for further details)	    		        
		    
		  */
			void setWeights(const std::vector<Int>& weight_labels, const std::vector<DoubleReal>& weights);

	 private:
     /**
        @brief find next grid search parameter combination

        The current grid cell is given in @p actual_values.
        The result is returned in @p actual_values.

     */
      bool nextGrid_(const std::vector<DoubleReal>& start_values, 
                     const std::vector<DoubleReal>& step_sizes,
                     const std::vector<DoubleReal>& end_values,
                     const bool additive_step_sizes,
                     std::vector<DoubleReal>& actual_values);

			Size getNumberOfEnclosedPoints_(DoubleReal m1, DoubleReal m2, const std::vector<std::pair<DoubleReal, DoubleReal> >& 	points);
	
		  /**
		    @brief Initializes the svm with standard parameters
		    
		  */
	    void initParameters_();

      /**
        @brief This function is passed to lib svm for output control

        The intention is to discard the output, as we don't need it.

      */
      static void printToVoid(const char * /*s*/);

	    svm_parameter* 												param_;  	       	    // the parameters for the svm
	    svm_model*     												model_;   			      // the learnt svm discriminant
	    DoubleReal 														sigma_;								// for the oligo kernel (amount of positional smearing) 
			std::vector<DoubleReal>								sigmas_;							// for the combined oligo kernel (amount of positional smearing) 
			std::vector<DoubleReal>								gauss_table_;					// lookup table for fast computation of the oligo kernel
			std::vector<std::vector<DoubleReal>	> gauss_tables_;				// lookup table for fast computation of the combined oligo kernel
			Size			 														kernel_type_;					// the actual kernel type	
			Size			 														border_length_;				// the actual kernel type				
			svm_problem*													training_set_;				// the training set
			svm_problem*													training_problem_;		// the training set
			SVMData 															training_data_;				// the training set (different encoding)

	};
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_SVM_SVMWRAPPER_H