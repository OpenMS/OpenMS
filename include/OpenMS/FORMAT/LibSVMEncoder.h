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

#ifndef OPENMS_FORMAT_LIBSVMENCODER_H
#define OPENMS_FORMAT_LIBSVMENCODER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <svm.h>

#include <vector>
#include <utility>

namespace OpenMS 
{
  /**
    @brief Serves for encoding sequences into feature vectors
    
    The class can be used to construct composition vectors for
    sequences. Additionally the vectors can be encoded into
    the libsvm format.
    
  */
  class LibSVMEncoder
  {
    public:
      /// Constructor
      LibSVMEncoder();
      /// Destructor
      ~LibSVMEncoder();
            
      /**
 				@brief returns a composition vector of 'sequence'
 				
 				The allowed characters given by 'allowed_characters' are counted in the sequence 'sequence'
 				and the relative frequency of the letters are sored in the composition vector. 
 				The first entry of the vector (<UnsignedInt, DoubleReal>) corresponds to the first letter of 
 				'allowed_characters' that has a non zero frequency in 'sequence' and its corresponding 
 				relative frequency...
			*/ 				
      std::vector< std::pair<SignedInt, DoubleReal> >* encodeCompositionVector(const String& sequence, 
																												 								 			   const String& allowed_characters = "ACDEFGHIKLMNPQRSTVWY");      
																																		  
      /**
 				@brief returns composition vector of the sequences given by 'sequence'
 				
 				The allowed characters given by 'allowed_characters' are counted in the sequences 'sequences'
 				and the relative frequency of the letters are sored in the composition vectors. 
 				The first entry of the first vector (<UnsignedInt, DoubleReal>) corresponds to the first letter of 
 				'allowed_characters' that has a non zero frequency in the first 'sequence' and its corresponding 
 				relative frequency...
			*/ 				
      std::vector< std::vector< std::pair<SignedInt, DoubleReal> > >* encodeCompositionVectors(const std::vector<String>& sequences, 
																												 								 			   						const String& allowed_characters);      
			/// encodes the feature vector in LibSVM compliant format																																		  
      svm_node* encodeLibSVMVector(
      	const std::vector< std::pair<SignedInt, DoubleReal> >& feature_vector);
      
			/// encodes the feature vectors in LibSVM compliant format																																		  
      std::vector<svm_node*>* encodeLibSVMVectors(
      	const std::vector< std::vector< std::pair<SignedInt, DoubleReal> > >& feature_vectors);
      
			/// encodes the LibSVM compliant vectors into a LibSVM compliant structure																																		  
      svm_problem* encodeLibSVMProblem(const std::vector<svm_node*>&  vectors, 
      																 std::vector<DoubleReal>* 		  labels);
      
      /// creates composition vectors for 'sequences' and stores them in LibSVM compliant format
			svm_problem* encodeLibSVMProblemWithCompositionVectors(const std::vector<String>& sequences,
																														 std::vector<DoubleReal>*   labels,
																														 const String&              allowed_characters);      
    
      /// creates composition vectors with additional length information for 'sequences' and stores them in LibSVM compliant format
			svm_problem* encodeLibSVMProblemWithCompositionAndLengthVectors(const std::vector<String>& sequences,
																																			std::vector<DoubleReal>*   labels,
																																			const String&              allowed_characters,
																																			UnsignedInt                maximum_sequence_length);      
    	
    	/// stores the LibSVM-encoded data in a text file that can be used by the LibSVM applications (svm-scale, svm-train,...)
			bool storeLibSVMProblem(const String& filename, const svm_problem* problem) const;

    	/// loads the LibSVM-encoded data stored in 'filename'
			svm_problem* loadLibSVMProblem(const String& filename);

    	/// encodes the borders of the sequence as k_mer oligos and stores them in 'libsvm_vector'
			void encodeOligoBorders(String 																						 sequence,
															UnsignedInt 																			 k_mer_length,
															const String& 																		 allowed_characters,
															UnsignedInt                                        border_length,
															std::vector< std::pair<SignedInt, DoubleReal> >& libsvm_vector,
															bool 																							 strict = false,
															bool 																							 length_encoding = false);

      /// creates oligo border vectors vectors for 'sequences' and stores them in LibSVM compliant format
			svm_problem* encodeLibSVMProblemWithOligoBorderVectors(const std::vector<String>&     sequences,
																														 std::vector<DoubleReal>*  			labels,
																														 UnsignedInt 										k_mer_length,
																														 const String& 	 				  			allowed_characters,
																														 UnsignedInt 										border_length,
																											 			 bool 													strict = false,
																											 			 bool 													length_encoding = false);

      /// creates combined oligo border vectors vectors for 'sequences' and stores them in LibSVM compliant format
			svm_problem* encodeLibSVMProblemWithCombinedOligoBorderVectors(const std::vector<String>&     sequences,
																														 std::vector<DoubleReal>*  			labels,
																														 const std::vector<std::pair<UnsignedInt, UnsignedInt> >& parameters,
																														 const std::vector<DoubleReal>& sigmas,
																														 const String& 	 				  			allowed_characters,
																														 bool 													strict = false,
																											 			 bool 													length_encoding = false);

			/// stores a string representation of 'vector' in 'output'
			void oligoBorderVectorToString(svm_node* vector, UnsignedInt border_length, String& output);

      /**
 				@brief stores a string representation of the encoded sequences in 'vectors' in 'output'
 				
 				This function can be used if one wants to print the feature vectors that are used in
 				the oligo kernel. 'border_length' is the parameter with which the sequence had been encoded.
			*/ 				
			void oligoBorderVectorsToString(svm_problem* vector, UnsignedInt border_length, String& output);
																																	 
      /**
 				@brief stores a string representation of the encoded sequence 'vector' in 'output'
 				
 				This function can be used if one wants to print one feature vector that is used in
 				the libsvm.
			*/ 				
			void libSVMVectorToString(svm_node* vector, String& output);

      /**
 				@brief stores a string representation of the encoded sequences in 'vectors' in 'output'
 				
 				This function can be used if one wants to print the feature vectors that are used in
 				the libsvm.
			*/ 				
			void libSVMVectorsToString(svm_problem* vector, String& output);

			typedef DoubleReal (*double_pt_2_string_double)(const String&, DoubleReal);				

			static DoubleReal getPeptideCharge(const String& sequence, DoubleReal ph);
				
			static DoubleReal getPeptideWeight(const String& sequence, DoubleReal charge = 0);

			static DoubleReal getPeptideSequenceIndex(const String& sequence, DoubleReal scale = 1);
		
			void encodeVector(const String& sequence, DoubleReal parameter, std::vector<double_pt_2_string_double> functions, std::vector< std::pair<SignedInt, DoubleReal> >& encoded_vector, UnsignedInt start_index = 1);	
			

  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_LIBSVMENCODER_H
