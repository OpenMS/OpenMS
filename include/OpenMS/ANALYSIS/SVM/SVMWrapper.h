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

#ifndef OPENMS_ANALYSIS_SVM_SVMWRAPPER_H
#define OPENMS_ANALYSIS_SVM_SVMWRAPPER_H

#include <svm.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <string>
#include <vector>
#include <map>
#include <math.h>

namespace OpenMS 
{
	
  /**
    @brief 	Parameters for the svm to be set from outside

		This type is used to specify the kind of parameter that
		is to be set or retrieved by the set/getParameter methods.
		
		SVM_TYPE:  		the svm type cab be NU_SVR or EPSILON_SVR
	 
	  KERNEL_TYPE: 	the kernel type
	 
   	DEGREE:      	the degree for the polynomial- kernel
	 
   	C:           	the C parameter of the svm
   	NU:          	the nu parameter for nu-SVR
   	P:           	the epsilon parameter for epsilon-SVR
   	GAMMA:       	the gamma parameter of the POLY, RBF and SIGMOID kernel
	*/
	typedef enum{
		
		SVM_TYPE,
	  KERNEL_TYPE,
	  DEGREE,
	  C,
	  NU,
	  P,
	  GAMMA,
	  PROBABILITY,
	  SIGMA,
	  BORDER_LENGTH,
	  SIGMA_1,
	  SIGMA_2,
	  SIGMA_3,
	  SIGMA_4,
	  BORDER_LENGTH_1,
	  BORDER_LENGTH_2,
	  BORDER_LENGTH_4,
	  BORDER_LENGTH_3
	  	    	
	}SVM_parameter_type;
	
	typedef enum 
	{ 
		OLIGO = 19,
		OLIGO_COMBINED
	}SVM_kernel_type; /* kernel_type */
	
	class SVMWrapper
	{
	 public:
	
			/// standard constructor
	    SVMWrapper();

			/// destructor	
	    ~SVMWrapper();
	
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
	    void setParameter(SVM_parameter_type type, int value);
	
		  /**
		    @brief sets the double parameters of the svm
		    
			*/
	    void setParameter(SVM_parameter_type type, double value);
		
		  /**
		    @brief	trains the svm 

	      The svm is trained with the data stored in the 'svm_problem' structure.
			*/
	    int  train(struct svm_problem* problem);
	
		  /**
		    @brief	saves the svm model 

	      The model of the trained svm is saved into 'modelFilename'.
			*/
	    void saveModel(std::string modelFilename);
	
		  /**
		    @brief loads the model

	      The svm- model is loaded. After this, the svm is ready for 
	      prediction.
		  */
	    void loadModel(std::string modelFilename);
	
		  /**
		    @brief predicts the labels using the trained model
		    
	     	 The prediction process is started and the results are returned.

		  */
	    std::vector<DoubleReal>* predict(struct svm_problem* predictProblem);

		  /**
		    @brief predicts the labels using the trained model
		    
	     	 The prediction process is started and the results are returned.

		  */
	    std::vector<DoubleReal>* predict(const std::vector<svm_node*>& vectors);
	
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
	    int getIntParameter(SVM_parameter_type type);
	
		  /**
		    @brief You can get the actual double- parameters of the svm
		    
		    C:            the C parameter of the svm
		    P:			      the P parameter of the svm (sets the epsilon in
	     											  epsilon-svr)
	   		NU:           the nu parameter in nu-SVR
	   		GAMMA:				for POLY, RBF and SIGMOID		    
		  */
	    double getDoubleParameter(SVM_parameter_type type); 

		  /**
		    @brief You can create 'number' equally sized random partitions 
		    
		  */
			static std::vector<svm_problem*>* createRandomPartitions(svm_problem* problem,
															 											 		UnsignedInt  number);
	
		  /**
		    @brief You can merge partitions excuding the partition with index 'except' 
		    
		  */
	    static svm_problem* mergePartitions(const std::vector<svm_problem*>* const problems,
																	 				UnsignedInt 								except);
																	 				
		  /**
		    @brief Returns the stored labels of the encoded SVM data 
		    
		  */
			static std::vector<DoubleReal>* getLabels(svm_problem* problem);
																	 				
		  /**
		    @brief Returns the stored labels of the encoded SVM data 
		    
		  */
			std::map<SVM_parameter_type, DoubleReal>* performCrossValidation(svm_problem* problem,
																					std::map<SVM_parameter_type, DoubleReal>& start_values,
																					std::map<SVM_parameter_type, DoubleReal>& step_sizes,
																					std::map<SVM_parameter_type, DoubleReal>& end_values,
																 					DoubleReal* cv_quality,
																 					UnsignedInt number_of_partitions,
																 					UnsignedInt number_of_runs,
																 					bool 			  additive_step_size = true,
																 					bool output = false,
																 					String 		  performances_file_name = "performances.txt");
																 					
		  /**
		    @brief Returns the probability parameter sigma of the fitted laplace model.		      
		    
		    The libsvm is used to fit a laplace model to the prediction values by performing
		    an internal cv using the training set if setParameter(PROBABILITY, 1) was invoked
		    before using train. Look for your libsvm documentation for more details.
		    The model parameter sigma is returned by this method.	If no model was fitted during 
		    training zero is returned. 
		  */
			double getSVRProbability();																			 					

		  /**
		    @brief calculates the oligo kernel value for the encoded sequences 'x' and 'y'
		    
		    This kernel function calculates the oligo kernel value [Meinicke 04] for
		    the sequences 'x' and 'y' that had been encoded by the encodeOligoBorder... function
		    of the LibSVMEncoder class. 	      		    
		  */
			static DoubleReal kernelOligo(const svm_node*									x, 
																	  const svm_node*									y,
																	  const std::vector<DoubleReal>&	gauss_table,																	  
								  								  DoubleReal                      sigma_square = 0,	
								  								  UnsignedInt	 										max_distance = 50);

		  /**
		    @brief calculates the significance borders of the error model and stores them in 'borders'	      
		    
		  */
			void getSignificanceBorders(svm_problem* data, 
																	std::pair<DoubleReal, DoubleReal>& borders,
																	DoubleReal confidence = 0.95,
																	UnsignedInt number_of_runs = 10,
																	UnsignedInt number_of_partitions = 5,
																	DoubleReal step_size = 0.01,
																	UnsignedInt max_iterations = 1000000);

		  /**
		    @brief calculates a p-value for a given data point using the model parameters
		    
		    Uses the model parameters to calculate the p-value for 'point' which has the data
		    entries: measured, predicted retention time.	      
		    
		  */
			DoubleReal getPValue(DoubleReal 													sigma1, 
													 DoubleReal 													sigma2,
													 std::pair<DoubleReal, DoubleReal>		point);

		  /**
		    @brief stores the prediction values for the encoded data in 'decision_values'
		    
		    This function can be used to get the prediction values of the data if a model 
		    is already trained by the train() method.
		    
		  */
			void getDecisionValues(svm_problem* data, std::vector<DoubleReal>& decision_values);

		  /**
		    @brief Scales the data such that every coloumn is scaled to [-1, 1].
		    
		  */
			void scaleData(svm_problem* data, UnsignedInt number_of_combinations = 1, UnsignedInt start_combination = 0, SignedInt max_scale_value = -1);

  		void setTrainingSample(svm_problem* training_sample);

			static void calculateGaussTable(UnsignedInt 							border_length, 
																			 DoubleReal 							sigma, 
																			 std::vector<DoubleReal>&	gauss_table);

		  /**
		    @brief computes the kernel matrix using the actual svm parameters	and the given data
		    
		    This function can be used to compute a kernel matrix. 'problem1' and 'problem2'
		    are used together wit the oligo kernel function (could be extended if you 
		    want to use your own kernel functions).      
		    
		  */
			svm_problem* computeKernelMatrix(svm_problem* problem1, svm_problem* problem2);

		protected:
						
			void destroyProblem(svm_problem* problem);
																	 				
	 private:
			UnsignedInt getNumberOfEnclosedPoints(DoubleReal m1, 
																						DoubleReal m2, 
																						const std::vector<std::pair<DoubleReal, DoubleReal> >& 	points);
	
	    svm_parameter* 												param_;  	       	    // the parameters for the svm
	    svm_model*     												model_;   			      // the learnt svm discriminant
	    DoubleReal 														sigma_;								// for the oligo kernel (amount of positional smearing) 
			std::vector<DoubleReal>								sigmas_;							// for the combined oligo kernel (amount of positional smearing) 
			std::vector<DoubleReal>								gauss_table_;					// lookup table for fast computation of the oligo kernel
			std::vector<std::vector<DoubleReal>	> gauss_tables_;				// lookup table for fast computation of the combined oligo kernel
			UnsignedInt			 											kernel_type_;					// the actual kernel type	
			UnsignedInt			 											border_length_;				// the actual kernel type				
			svm_problem*													training_set_;				// the training set
			svm_problem*													training_problem_;		// the training set
			bool 																	loaded_model_;				// indicates wether a model is loaded from a file or not

		  /**
		    @brief Initializes the svm with standard parameters
		    
		  */
	    void initParameters();
	};
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_SVM_SVMWRAPPER_H
