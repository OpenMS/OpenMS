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
      /// Copy constructor
      LibSVMEncoder(const LibSVMEncoder& source);
      /// Destructor
      ~LibSVMEncoder();
      
      /// Assignment operator
      LibSVMEncoder& operator = (const LibSVMEncoder& source);
      
      std::vector< std::pair<UnsignedInt, DoubleReal> >* encodeCompositionVector(const String& sequence, 
																												 								 			   const String& allowed_characters);      
																																		  
      std::vector< std::vector< std::pair<UnsignedInt, DoubleReal> > >* encodeCompositionVectors(const std::vector<String>& sequences, 
																												 								 			   						const String& allowed_characters);      
																																		  
      svm_node* encodeLIBSVMVector(
      	const std::vector< std::pair<UnsignedInt, DoubleReal> >& feature_vector);
      
      std::vector<svm_node*>* encodeLIBSVMVectors(
      	const std::vector< std::vector< std::pair<UnsignedInt, DoubleReal> > >& feature_vectors);
      
      svm_problem* encodeLIBSVMProblem(const std::vector<svm_node*>&  vectors, 
      																 std::vector<DoubleReal>* 		  labels);
      
			svm_problem* encodeLIBSVMProblemWithCompositionVectors(const std::vector<String>& sequences,
																														 std::vector<DoubleReal>*   labels,
																														 const String&              allowed_characters);      
    
			svm_problem* encodeLIBSVMProblemWithCompositionAndLengthVectors(const std::vector<String>& sequences,
																																			std::vector<DoubleReal>*   labels,
																																			const String&              allowed_characters,
																																			UnsignedInt                maximum_sequence_length);      
    
    protected:


  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_LIBSVMENCODER_H
