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

#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <map>
#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS 
{
  LibSVMEncoder::LibSVMEncoder() 
  {
    
  }
  
  LibSVMEncoder::~LibSVMEncoder()
  {
    
  }
		
	vector< pair<SignedInt, DoubleReal> >* LibSVMEncoder::encodeCompositionVector(const String& sequence, 
																											 				 					    const String& allowed_characters)
	{
	
		UnsignedInt number_of_different_letters = allowed_characters.size();
		UnsignedInt* counts = new UnsignedInt[number_of_different_letters];
		UnsignedInt total_count = 0;
		vector< pair<SignedInt, DoubleReal> >* composition_vector = 
			new vector< pair<SignedInt, DoubleReal> >();
		
		for (UnsignedInt i = 0; i < number_of_different_letters; i++)
		{
			counts[i] = 0;
		}
		
		for (UnsignedInt i = 0; i < sequence.size(); i++)
		{			
			if (allowed_characters.find(sequence[i]) != String::npos)
			{			
				counts[allowed_characters.find(sequence[i])]++;
				total_count++;
			}
		}
		
		for(UnsignedInt i = 0; i < number_of_different_letters; i++)
		{
			if (counts[i] > 0)
			{
				composition_vector->push_back(make_pair(i + 1, (((DoubleReal) counts[i]) / total_count)));
			}
		}
		delete [] counts;
		
		return composition_vector;
	}
		
	vector< vector< pair<SignedInt, DoubleReal> > >* LibSVMEncoder::encodeCompositionVectors(const vector<String>& sequences, 
																												 								 			   						 const String& allowed_characters)
	{
		vector< vector< pair<SignedInt, DoubleReal> > >* composition_vectors = 
			new vector< vector< pair<SignedInt, DoubleReal> > >();
		vector< pair<SignedInt, DoubleReal> >* composition_vector;
		
		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{
			composition_vector = encodeCompositionVector(sequences[i],
																									 allowed_characters);
			composition_vectors->push_back(*composition_vector);
			delete composition_vector;
		}
					
		return composition_vectors;
	}      
	
	svm_node* LibSVMEncoder::encodeLibSVMVector(const vector< pair<SignedInt, DoubleReal> >& feature_vector)
	{
		
		vector< pair<SignedInt, DoubleReal> >::const_iterator vector_iterator;
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
	
	vector<svm_node*>* LibSVMEncoder::encodeLibSVMVectors(
  	const vector< vector< pair<SignedInt, DoubleReal> > >& feature_vectors)
  {
  	vector<svm_node*>* libsvm_vectors = new vector<svm_node*>();
  	
  	for(UnsignedInt i = 0; i < feature_vectors.size(); i++)
  	{
  		libsvm_vectors->push_back(encodeLibSVMVector(feature_vectors[i]));
  	}
  	
  	return libsvm_vectors;
  }
	
	
	svm_problem* LibSVMEncoder::encodeLibSVMProblem(const vector<svm_node*>& vectors, 
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
		
//		cout << "Problem encoded" << endl;
		
		return problem; 
	}
	
	svm_problem* LibSVMEncoder::encodeLibSVMProblemWithCompositionVectors(const vector<String>&    sequences,
																																	std::vector<DoubleReal>* labels,
																																	const String&   				 allowed_characters)
	{
		vector<svm_node*> vectors;
		vector< pair<SignedInt, DoubleReal> >* encoded_vector;
		svm_node* libsvm_vector;
		
		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{
			
			encoded_vector = encodeCompositionVector(sequences[i], allowed_characters);
			libsvm_vector = encodeLibSVMVector(*encoded_vector);
			vectors.push_back(libsvm_vector);
		}
		
		return encodeLibSVMProblem(vectors, labels);		
	} 

	svm_problem* LibSVMEncoder::encodeLibSVMProblemWithCompositionAndLengthVectors(const vector<String>&    sequences,
																																								 std::vector<DoubleReal>* labels,
																																								 const String& 	 				  allowed_characters,
																																								 UnsignedInt 							maximum_sequence_length)
	{
		vector<svm_node*> vectors;
		vector< pair<SignedInt, DoubleReal> >* encoded_vector;
		svm_node* libsvm_vector;
		
		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{
			
			encoded_vector = encodeCompositionVector(sequences[i], allowed_characters);
			encoded_vector->push_back(make_pair(allowed_characters.size() + 1, ((DoubleReal) sequences[i].length()) / maximum_sequence_length));
			libsvm_vector = encodeLibSVMVector(*encoded_vector);
			vectors.push_back(libsvm_vector);
		}
		
		return encodeLibSVMProblem(vectors, labels);		
	}
	
	bool LibSVMEncoder::storeLibSVMProblem(const String& filename, const svm_problem* problem) const
	{
		ofstream output_file(filename.c_str());
		SignedInt j = 0;
		SignedInt counter = 0;
		
		if (problem == NULL)
		{
			return false;
		}
		
		// checking if file is writable
		if (!File::writable(filename))
		{
			return false;
		}
			
		// writing feature vectors		
		for(SignedInt i = 0; i < problem->l; i++)
		{
			j = 0;
			counter = 0;
			output_file << problem->y[i] << " ";
			while(problem->x[i][j].index != -1)
			{
				output_file << problem->x[i][j].index << ":" << problem->x[i][j].value << " " ;					
				++j;
			}
			output_file << endl;
		}
		output_file.flush();
		output_file.close();
		cout.flush();
		return true;
	}
	
	svm_problem* LibSVMEncoder::loadLibSVMProblem(const String& filename)
	{
		svm_problem* data = NULL;				
		UnsignedInt counter = 0;
		vector<String> parts;
		vector<String> temp_parts;
		
		if (!File::exists(filename))
		{
			return NULL;
		}
		if (!File::readable(filename))
		{
			return NULL;
		}
    if (File::empty(filename))
    {
			return NULL;
    }		

		TextFile text_file(filename.c_str(), true);
    TextFile::iterator it;

		it = text_file.begin();
		
		data = new svm_problem;
		data->l = text_file.size();
		data->x = new svm_node*[text_file.size()];
		data->y = new DoubleReal[text_file.size()];
		while(counter < text_file.size()&& it != text_file.end())
		{
			it->split(' ', parts);
			data->y[counter] = parts[0].trim().toFloat();		
			data->x[counter] = new svm_node[parts.size()];			
			for(UnsignedInt j = 1; j < parts.size(); ++j)
			{
				parts[j].split(':', temp_parts);
				if (temp_parts.size() < 2)
				{
					delete data;
					return NULL;
				}
				data->x[counter][j - 1].index = temp_parts[0].trim().toInt();
				data->x[counter][j - 1].value = temp_parts[1].trim().toFloat();
			}
			data->x[counter][parts.size() - 1].index = -1;
			data->x[counter][parts.size() - 1].value = 0;
			++counter;
			++it;
		}
		return data;				
	}

	void LibSVMEncoder::encodeOligoBorders(String sequence,
													UnsignedInt k_mer_length,
													const String& allowed_characters,
													UnsignedInt border_length,
													vector< pair<SignedInt, DoubleReal> >& libsvm_vector,
													bool strict,
													bool length_encoding)
	{
	  multimap<SignedInt, SignedInt>  	         			ordered_tree;
	  multimap<SignedInt, SignedInt>::iterator 				elements;
	  pair<SignedInt, SignedInt>  	             			values;
	  UnsignedInt               										   	oligo_value = 0;
	  UnsignedInt                										   	factor      = 1;
	  map<String::value_type, UnsignedInt>							residue_values;
	  UnsignedInt 																			counter 		= 0;
	  UnsignedInt 																			number_of_residues = 0;
	  UnsignedInt 																			left_border = 0;
	  UnsignedInt																				right_border = 0;
	  UnsignedInt																				sequence_length = 0;

	  number_of_residues = allowed_characters.size();
	  
  	libsvm_vector.clear();
	  sequence_length = sequence.size();
	  if (k_mer_length <= sequence_length)
	  {
	  	// if a border must not be longer than half of the peptide
	  	if (strict)
	  	{
			  if (border_length > (sequence_length - k_mer_length + 1) / 2)
			  {
			  	left_border = (UnsignedInt) (floor((sequence_length - k_mer_length + 1) / 2));
			  	right_border = (UnsignedInt) (ceil((sequence_length - k_mer_length + 1) / 2));
			  }
			  else
			  {
			  	left_border = border_length;
			  	right_border = sequence_length - k_mer_length + 1 - border_length;
			  }
	  	}
	  	else
	  	{
			  if (border_length >= sequence_length - k_mer_length + 1)
			  {
			  	left_border = sequence_length - k_mer_length + 1;
			  	right_border = 0;
			  }
			  else
			  {
			  	left_border = border_length;
			  	right_border = sequence_length - k_mer_length + 1 - border_length;
			  }
			}	

			for(UnsignedInt i = 0; i < number_of_residues; ++i)
			{	
				residue_values.insert(make_pair(allowed_characters[i], counter));
				++counter;
			}
		  for(int k = k_mer_length - 1; k >= 0; k--)
		  {
				oligo_value += factor * residue_values[sequence[k]];
				factor *= number_of_residues;
			}
		  factor /= number_of_residues;
		  values.first = ((SignedInt) (oligo_value + 2));
		  values.second = 1;
		  ordered_tree.insert(values);
		
		  for(UnsignedInt j = 1; j < left_border; j++)
		  {
				oligo_value -= factor * residue_values[sequence[j - 1]];
				oligo_value = oligo_value * number_of_residues + residue_values[sequence[j + k_mer_length - 1]];
		
				values.first = ((SignedInt) (oligo_value + 2));
				values.second = j + 1;
		
				ordered_tree.insert(values);	
		  }
		  oligo_value = 0;
		  factor = 1;
		  
		  if (k_mer_length > 1)
		  {
			  for(int k = k_mer_length; k > 0; k--)
			  {
					oligo_value += factor * residue_values[sequence[sequence_length - k]];
					factor *= number_of_residues;
				}
			  factor /= number_of_residues;
			  values.first = ((SignedInt) (oligo_value + 2));
			  values.second = 1;
			  ordered_tree.insert(values);
			
			  for(UnsignedInt j = 1; j < left_border; j++)
			  {
					oligo_value -= factor * residue_values[sequence[sequence_length - j]];
					oligo_value = oligo_value * number_of_residues + residue_values[sequence[sequence_length - k_mer_length - j]];
			
					values.first = ((SignedInt) (oligo_value + 2));
					values.second = j + 1;
			
					ordered_tree.insert(values);	
			  }
		  }
		  else
		  {
			  for(UnsignedInt k = right_border + k_mer_length; k > right_border; k--)
			  {
					oligo_value += factor * residue_values[sequence[k - 1]];
					factor *= number_of_residues;
				}
			  factor /= number_of_residues;
			  values.first = oligo_value + 2;
			  values.second = (right_border - sequence_length) * -1;
			  ordered_tree.insert(values);
			  for(UnsignedInt j = right_border + 1; j < sequence_length - k_mer_length + 1; j++)
			  {
					oligo_value -= factor * residue_values[sequence[j - 1]];
					oligo_value = oligo_value * number_of_residues + residue_values[sequence[j + k_mer_length - 1]];
			
					values.first = oligo_value + 2;
					values.second = (j - sequence_length) * -1;
			
					ordered_tree.insert(values);	
			  }
			}		  	
		  for(elements = ordered_tree.begin(); elements != ordered_tree.end(); ++elements)
		  {
				libsvm_vector.push_back(make_pair(elements->first, elements->second));	
		  }
		  if (length_encoding)
		  {
		  	libsvm_vector.push_back(make_pair((SignedInt) sequence.size(), pow(k_mer_length, number_of_residues) + 1));
		  }
		}

	}
	
	svm_problem* LibSVMEncoder::encodeLibSVMProblemWithOligoBorderVectors(const vector<String>&     sequences,
																																				vector<DoubleReal>*       labels,
																																				UnsignedInt 							k_mer_length,
																																				const String& 	 				  allowed_characters,
																																				UnsignedInt 							border_length,
																																				bool											strict,
																																				bool 											length_encoding)
	{
		vector<svm_node*> vectors;
		vector< pair<SignedInt, DoubleReal> > encoded_vector;
		svm_node* libsvm_vector;
		
		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{			
			encodeOligoBorders(sequences[i], k_mer_length, allowed_characters, border_length, encoded_vector, strict, length_encoding);
			libsvm_vector = encodeLibSVMVector(encoded_vector);
			vectors.push_back(libsvm_vector);
		}
		
		return encodeLibSVMProblem(vectors, labels);		
	}
	
	void LibSVMEncoder::libSVMVectorToString(svm_node* vector, String& output)
	{
		UnsignedInt i = 0;
		
		output.clear();
		while (vector[i].index != -1)
		{
			output = output + "(" + String(vector[i].index) + ", " + String(vector[i].value) + ") ";
			++i;
		}
	}
	
	void LibSVMEncoder::libSVMVectorsToString(svm_problem* vector, String& output)
	{
		String temp_string = "";
		
		output.clear();
		if (vector != NULL)
		{
			for(SignedInt i = 0; i < vector->l; ++i)
			{
				libSVMVectorToString(vector->x[i], temp_string);
				output = output + temp_string + "\n";
				temp_string = "";
			}
		}
	}
	
} // namespace OpenMS
