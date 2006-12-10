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

#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <iostream>

using namespace std;

namespace OpenMS 
{
  LibSVMEncoder::LibSVMEncoder() 
  {
    
  }
  
  LibSVMEncoder::~LibSVMEncoder()
  {
    
  }
		
	vector< pair<UnsignedInt, DoubleReal> >* LibSVMEncoder::encodeCompositionVector(const String& sequence, 
																											 				 					    const String& allowed_characters)
	{
	
		UnsignedInt number_of_different_letters = allowed_characters.size();
		UnsignedInt* counts = new UnsignedInt[number_of_different_letters];
		UnsignedInt total_count = 0;
		vector< pair<UnsignedInt, DoubleReal> >* composition_vector = 
			new vector< pair<UnsignedInt, DoubleReal> >();
		
		for (UnsignedInt i = 0; i < number_of_different_letters; i++)
		{
			counts[i] = 0;
		}
		
		for (UnsignedInt i = 0; i < sequence.size(); i++)
		{			
			if (allowed_characters.find(sequence[i]) != string::npos)
			{			
				counts[allowed_characters.find(sequence[i])]++;
				total_count++;
			}
		}
		
		for(UnsignedInt i = 0; i < number_of_different_letters; i++)
		{
			if (counts[i] > 0)
			{
				composition_vector->push_back(make_pair(i, (((DoubleReal) counts[i]) / total_count)));
			}
		}
		delete [] counts;
		
		return composition_vector;
	}
	
	vector< vector< pair<UnsignedInt, DoubleReal> > >* LibSVMEncoder::encodeCompositionVectors(const vector<String>& sequences, 
																												 								 			   		const String& allowed_characters)
	{
		vector< vector< pair<UnsignedInt, DoubleReal> > >* composition_vectors = 
			new vector< vector< pair<UnsignedInt, DoubleReal> > >();
		vector< pair<UnsignedInt, DoubleReal> >* composition_vector;
		
		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{
			composition_vector = encodeCompositionVector(sequences[i],
																									 allowed_characters);
			composition_vectors->push_back(*composition_vector);
			delete composition_vector;
		}
					
		return composition_vectors;
	}      
	
	
	svm_node* LibSVMEncoder::encodeLIBSVMVector(const vector< pair<UnsignedInt, DoubleReal> >& feature_vector)
	{
		
		vector< pair<UnsignedInt, DoubleReal> >::const_iterator vector_iterator;
		svm_node* nodes;
		UnsignedInt i = 0;

		nodes = new svm_node[feature_vector.size() + 1];					   
		
		for(vector_iterator = feature_vector.begin(); 
				vector_iterator != feature_vector.end(); 
				vector_iterator++)
		{
			nodes[i].index = vector_iterator->first;
			nodes[i].value = vector_iterator->second;
			i++;
		}
		nodes[feature_vector.size()].index = -1;
		nodes[feature_vector.size()].value = 0;
		
		return nodes;	
	}
	
	vector<svm_node*>* LibSVMEncoder::encodeLIBSVMVectors(
  	const vector< vector< pair<UnsignedInt, DoubleReal> > >& feature_vectors)
  {
  	vector<svm_node*>* libsvm_vectors = new vector<svm_node*>();
  	
  	for(UnsignedInt i = 0; i < feature_vectors.size(); i++)
  	{
  		libsvm_vectors->push_back(encodeLIBSVMVector(feature_vectors[i]));
  	}
  	
  	return libsvm_vectors;
  }
	
	
	svm_problem* LibSVMEncoder::encodeLIBSVMProblem(const vector<svm_node*>& vectors, 
																						vector<DoubleReal>*  		 labels)
	{
		svm_problem* problem;
		svm_node** node_vectors;
		
		problem = new svm_problem;
		if (problem == NULL)
		{
			return NULL;
		}
		problem->l = (int) vectors.size();
		if (problem->l < 0)
		{
			return NULL;
		}
		
		problem->y = &((*labels)[0]);
		
		node_vectors = new svm_node*[problem->l];
		if (node_vectors == NULL)
		{
			free(problem);
			return NULL;
		}
		
		for(UnsignedInt i = 0; i < vectors.size(); i++)
		{
			node_vectors[i] = vectors[i];
		}
		problem->x = node_vectors;
		
		return problem; 
	}
	
	svm_problem* LibSVMEncoder::encodeLIBSVMProblemWithCompositionVectors(const vector<String>&    sequences,
																																	std::vector<DoubleReal>* labels,
																																	const String&   				 allowed_characters)
	{
		vector<svm_node*> vectors;
		vector< pair<UnsignedInt, DoubleReal> >* encoded_vector;
		svm_node* libsvm_vector;
		
		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{
			
			encoded_vector = encodeCompositionVector(sequences[i], allowed_characters);
			libsvm_vector = encodeLIBSVMVector(*encoded_vector);
			vectors.push_back(libsvm_vector);
		}
		
		return encodeLIBSVMProblem(vectors, labels);		
	} 

	svm_problem* LibSVMEncoder::encodeLIBSVMProblemWithCompositionAndLengthVectors(const vector<String>&    sequences,
																																								 std::vector<DoubleReal>* labels,
																																								 const String& 	 				  allowed_characters,
																																								 UnsignedInt 							maximum_sequence_length)
	{
		vector<svm_node*> vectors;
		vector< pair<UnsignedInt, DoubleReal> >* encoded_vector;
		svm_node* libsvm_vector;
		
		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{
			
			encoded_vector = encodeCompositionVector(sequences[i], allowed_characters);
			encoded_vector->push_back(make_pair(allowed_characters.size(), ((DoubleReal) sequences[i].length()) / maximum_sequence_length));
			libsvm_vector = encodeLIBSVMVector(*encoded_vector);
			vectors.push_back(libsvm_vector);
		}
		
		return encodeLIBSVMProblem(vectors, labels);		
	} 

} // namespace OpenMS
