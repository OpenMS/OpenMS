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
      std::vector< std::pair<UnsignedInt, DoubleReal> >* encodeCompositionVector(const String& sequence, 
																												 								 			   const String& allowed_characters);      
																																		  
      /**
 				@brief returns composition vector of the sequences given by 'sequence'
 				
 				The allowed characters given by 'allowed_characters' are counted in the sequences 'sequences'
 				and the relative frequency of the letters are sored in the composition vectors. 
 				The first entry of the first vector (<UnsignedInt, DoubleReal>) corresponds to the first letter of 
 				'allowed_characters' that has a non zero frequency in the first 'sequence' and its corresponding 
 				relative frequency...
			*/ 				
      std::vector< std::vector< std::pair<UnsignedInt, DoubleReal> > >* encodeCompositionVectors(const std::vector<String>& sequences, 
																												 								 			   						const String& allowed_characters);      
			/// encodes the feature vector in LIBSVM compliant format																																		  
      svm_node* encodeLIBSVMVector(
      	const std::vector< std::pair<UnsignedInt, DoubleReal> >& feature_vector);
      
			/// encodes the feature vectors in LIBSVM compliant format																																		  
      std::vector<svm_node*>* encodeLIBSVMVectors(
      	const std::vector< std::vector< std::pair<UnsignedInt, DoubleReal> > >& feature_vectors);
      
			/// encodes the LIBSVM compliant vectors into a LIBSVM compliant structure																																		  
      svm_problem* encodeLIBSVMProblem(const std::vector<svm_node*>&  vectors, 
      																 std::vector<DoubleReal>* 		  labels);
      
      /// creates composition vectors for 'sequences' and stores them in LIBSVM compliant format
			svm_problem* encodeLIBSVMProblemWithCompositionVectors(const std::vector<String>& sequences,
																														 std::vector<DoubleReal>*   labels,
																														 const String&              allowed_characters);      
    
      /// creates composition vectors with additional length information for 'sequences' and stores them in LIBSVM compliant format
			svm_problem* encodeLIBSVMProblemWithCompositionAndLengthVectors(const std::vector<String>& sequences,
																																			std::vector<DoubleReal>*   labels,
																																			const String&              allowed_characters,
																																			UnsignedInt                maximum_sequence_length);      
    
    protected:


  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_LIBSVMENCODER_H
