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
	*/
	typedef enum{
		
	    SVM_TYPE,
	    KERNEL_TYPE,
	    DEGREE,
	    C,
	    NU,
	    P,
	    PROBABILITY
	    	
	}SVM_parameter_type;
	
	class SVMWrapper
	{
	 public:
	
			/// standard constructor
	    SVMWrapper();

			/// destructor	
	    ~SVMWrapper();
	
		  /**
		    @brief You can set the parameters of the svm: 
	     
	         KERNEL_TYPE: can be POLY     for the polynomial kernel
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
		    P:			  the P parameter of the svm (sets the epsilon in
	     											  epsilon-svr)
	   		NU:           the nu parameter in nu-SVR		    
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
																 					UnsignedInt number_of_runs);
																 					
		  /**
		    @brief Returns the probability parameter sigma of the model.		      
		    
		  */
			double getSVRProbability();																			 					
																	 				

	 private:
	
	    struct svm_parameter* param_;               /// the parameters for the svm
	    struct svm_model*     model_;               /// the learnt svm discriminant 
	
		  /**
		    @brief Initializes the svm with standard parameters
		    
		  */
	    void initParameters();
	};
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_SVM_SVMWRAPPER_H







