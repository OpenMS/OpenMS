// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Sandro Andreotti $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_LIBSVMENCODER_H
#define OPENMS_FORMAT_LIBSVMENCODER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
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
  class OPENMS_DLLAPI LibSVMEncoder
  {
    public:
      /// Constructor
      LibSVMEncoder();
      /// Destructor
      ~LibSVMEncoder();
            
      /**
 				@brief stores a composition vector of 'sequence' in 'encoded_vector'
 				
 				The allowed characters given by 'allowed_characters' are counted in the sequence 'sequence'
 				and the relative frequency of the letters are sored in the composition vector. 
 				The first entry of the vector (<UInt, DoubleReal>) corresponds to the first letter of 
 				'allowed_characters' that has a non zero frequency in 'sequence' and its corresponding 
 				relative frequency...
			*/ 				
      void encodeCompositionVector(const String& sequence, std::vector< std::pair<Int, DoubleReal> >& encoded_vector, const String& allowed_characters = "ACDEFGHIKLMNPQRSTVWY");      
																																		  
      /**
 				@brief stores composition vectors of the sequences given by 'sequence' in 'composition_vectors'
 				
 				The allowed characters given by 'allowed_characters' are counted in the sequences 'sequences'
 				and the relative frequency of the letters are sored in the composition vectors. 
 				The first entry of the first vector (<UInt, DoubleReal>) corresponds to the first letter of 
 				'allowed_characters' that has a non zero frequency in the first 'sequence' and its corresponding 
 				relative frequency...
			*/ 				
      void encodeCompositionVectors(const std::vector<String>& sequences, const String& allowed_characters, std::vector< std::vector< std::pair<Int, DoubleReal> > >& composition_vectors);      
			/// encodes the feature vector in LibSVM compliant format																																		  
      svm_node* encodeLibSVMVector(const std::vector< std::pair<Int, DoubleReal> >& feature_vector);
      
			/// encodes the feature vectors in LibSVM compliant format																																		  
      void encodeLibSVMVectors(const std::vector< std::vector< std::pair<Int, DoubleReal> > >& feature_vectors, std::vector<svm_node*>& libsvm_vectors);
      
			/// encodes the LibSVM compliant vectors into a LibSVM compliant structure																																		  
      svm_problem* encodeLibSVMProblem(const std::vector<svm_node*>&  vectors, 
      																 std::vector<DoubleReal>& labels);
      
      /// creates composition vectors for 'sequences' and stores them in LibSVM compliant format
			svm_problem* encodeLibSVMProblemWithCompositionVectors(const std::vector<String>& sequences,
																														 std::vector<DoubleReal>&   labels,
																														 const String&              allowed_characters);      
    
      /// creates composition vectors with additional length information for 'sequences' and stores them in LibSVM compliant format
			svm_problem* encodeLibSVMProblemWithCompositionAndLengthVectors(const std::vector<String>& sequences,
																																			std::vector<DoubleReal>&   labels,
																																			const String&              allowed_characters,
																																			UInt                maximum_sequence_length);      
    	
      /// creates composition vectors with additional length and average weight information for 'sequences' and stores them in LibSVM compliant format
			svm_problem* encodeLibSVMProblemWithCompositionLengthAndWeightVectors(const std::vector<String>&    sequences,
																																						std::vector<DoubleReal>& labels,
																																						const String& 	 				  allowed_characters);

    	/// stores the LibSVM-encoded data in a text file that can be used by the LibSVM applications (svm-scale, svm-train,...)
			bool storeLibSVMProblem(const String& filename, const svm_problem* problem) const;

    	/// loads the LibSVM-encoded data stored in 'filename'
			svm_problem* loadLibSVMProblem(const String& filename);

    	/// encodes the borders of the sequence as k_mer oligos and stores them in 'libsvm_vector'
			void encodeOligoBorders(String 																						 sequence,
															UInt 																			 k_mer_length,
															const String& 																		 allowed_characters,
															UInt                                        border_length,
															std::vector< std::pair<Int, DoubleReal> >& libsvm_vector,
															bool 																							 strict = false,
															bool																							 unpaired = false,
															bool 																							 length_encoding = false);

      /// creates oligo border vectors vectors for 'sequences' and stores them in LibSVM compliant format
			svm_problem* encodeLibSVMProblemWithOligoBorderVectors(const std::vector<String>&     sequences,
																														 std::vector<DoubleReal>&  			labels,
																														 UInt                           k_mer_length,
																														 const String& 	 				  			allowed_characters,
																														 UInt                           border_length,
																											 			 bool 													strict = false,
																											 			 bool 													unpaired = false,
																											 			 bool 													length_encoding = false);

      /// creates oligo border vectors vectors for 'sequences' and stores them in 'vectors'
			void encodeProblemWithOligoBorderVectors(const std::vector<AASequence>&                             sequences,
																							 UInt                                                       k_mer_length,
																							 const String&                                              allowed_characters,
																							 UInt                                                       border_length,
																							 std::vector< std::vector< std::pair<Int, DoubleReal> > >& 	vectors);

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

      /**
 				@brief encodes an AASequence instance in oligo encoding
 				
 				This function is used to get the oligo encoding for AASequence 
        'sequence'. If a residue is modified, it gets an extra oligo function.
			*/ 				
			void encodeOligo(const AASequence& sequence,
								      UInt k_mer_length,
								      const String& allowed_characters,
								      std::vector< std::pair<Int, DoubleReal> >& values,
								      bool is_right_border = false);
								      
      /**
 				@brief frees all the memory of the svm_problem instance
 				
 				This function is used to free all the memory used by 'problem'
			*/ 				
			static void destroyProblem(svm_problem* problem);								      

		private:
			/// comparator for oligos encoded by encodeOligo
			static bool cmpOligos_(std::pair<Int, DoubleReal> a, 
                             std::pair<Int, DoubleReal> b);
											 
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_LIBSVMENCODER_H
