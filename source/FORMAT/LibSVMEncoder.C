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
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <map>
#include <iostream>
#include <fstream>

#include <qfile.h>
#include <qfileinfo.h>

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
	
	DoubleReal LibSVMEncoder::getPeptideCharge(const String& sequence, DoubleReal ph)
	{
		ResidueDB residue_db;
		DoubleReal sum = 0;
		DoubleReal temp = 0;
		const Residue* temp_residue;		
		
		for(String::const_iterator it = sequence.begin(); it != sequence.end(); ++it)			
		{
			temp_residue = residue_db.getResidue(*it);
			if (temp_residue->getOneLetterCode() == "E"
					|| temp_residue->getOneLetterCode() == "D")
			{
				temp = temp_residue->getPka();
				temp = pow(10, ph - temp);
				sum += (temp / (1.0 + temp));
			}
			else if (temp_residue->getOneLetterCode() == "H"
							|| temp_residue->getOneLetterCode() == "R"
							|| temp_residue->getOneLetterCode() == "K")
			{
				temp = temp_residue->getPka();
				temp = pow(10, ph - temp);
				sum += (1.0 / (1 + temp));
			}
			
		}
		return sum;		
		
	}																							
								
	DoubleReal LibSVMEncoder::getPeptideWeight(const String& sequence, DoubleReal charge)
	{
		ResidueDB residue_db;
		DoubleReal sum = 0;
		const Residue* temp_residue = NULL;		
		String::const_iterator it = sequence.begin();
			
		sum += residue_db.getResidue(*it)->getAverageWeight(Residue::NTerminal);
		++it;			
		
		for(; it != sequence.end(); ++it)			
		{
			temp_residue = residue_db.getResidue(*it);
			sum += temp_residue->getAverageWeight(Residue::Internal);						
		}
		sum -= temp_residue->getAverageWeight(Residue::Internal);	
		sum += temp_residue->getAverageWeight(Residue::NTerminal);
		return (sum + round(charge));	
	}

	DoubleReal LibSVMEncoder::getPeptideSequenceIndex(const String& sequence, DoubleReal scale)
	{
		ResidueDB residue_db;
		DoubleReal sum = 0;
		const Residue* temp_residue = NULL;	
		DoubleReal pi1 = 0;
		DoubleReal pi2 = 0;	
		DoubleReal temp_sum = 0;
		String::const_iterator it = sequence.begin();
			
		scale = 1;
			
		if (sequence.size() <= 1)
		{
			return 0.0;
		}
			
		temp_residue = residue_db.getResidue(*it);
		
		pi1 = temp_residue->getPiValue();
		++it;
		temp_residue = residue_db.getResidue(*it);
		pi2 = temp_residue->getPiValue();
		++it;
		for( ; it != sequence.end(); ++it)
		{
			pi1 = pi2;
			temp_residue = residue_db.getResidue(*it);
			pi2 = temp_residue->getPiValue();
			temp_sum = pi1 + pi2;
			sum = sum + temp_sum * temp_sum;			
		}
		return sqrt(sum / (sequence.size() - 1)) * scale;	
	}

	svm_node* LibSVMEncoder::encodeOHVector(const String&		sequence, 
																		 			DoubleReal 			ph)
	{
		vector<double_pt_2_string_double> functions;
		vector< pair<SignedInt, DoubleReal> > encoded_vector;
		
		functions.push_back(LibSVMEncoder::getPeptideWeight );
		functions.push_back(LibSVMEncoder::getPeptideCharge );
			
		encodeVector(sequence, ph, functions, encoded_vector);
		encoded_vector.push_back(make_pair(encoded_vector.size() + 1, sequence.size()));

		return encodeLibSVMVector(encoded_vector);		
	}

	void LibSVMEncoder::encodeVector(const String& sequence, DoubleReal parameter, vector<double_pt_2_string_double> functions, vector< pair<SignedInt, DoubleReal> >& encoded_vector, UnsignedInt start_index)
	{
		
		for(UnsignedInt i = 0; i < functions.size(); ++i)
		{
			encoded_vector.push_back(make_pair((i + start_index), (functions[i](sequence, parameter))));
		}
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
	
	svm_problem* LibSVMEncoder::encodeLibSVMProblemWithOHVectors(const std::vector<String>&    	sequences,
																																vector<DoubleReal>* 		 			labels,
																																DoubleReal 										ph)
	{
		vector<svm_node*> encoded_sequences;
		
		for(UnsignedInt i = 0; i < sequences.size(); ++i)
		{
			encoded_sequences.push_back(encodeOHVector(sequences[i], ph));
		}
		return encodeLibSVMProblem(encoded_sequences, labels);
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
	
	bool LibSVMEncoder::storeLibSVMProblem(const String& filename, const svm_problem* problem, SignedInt number_of_combinations) const
	{
		ofstream output_file(filename.c_str());
		QFile file;
		SignedInt j = 0;
		SignedInt counter = 0;
		
		if (problem == NULL)
		{
			return false;
		}
		
		// checking if file is writable
		file.setName(filename.c_str());
		file.open( IO_WriteOnly );
		if (!file.isWritable())
		{
			file.close();				
			return false;
		}
		file.close();
				
		// writing feature vectors		
		for(SignedInt i = 0; i < problem->l; i++)
		{
			j = 0;
			counter = 0;
			output_file << problem->y[i] << " ";
			if (number_of_combinations != -1)
			{
				while(counter < number_of_combinations * 2 - 1 || problem->x[i][j].index != -1)
				{
					if (problem->x[i][j].index == -1)
					{
						++counter;
					}
					if (number_of_combinations != 0 || problem->x[i][j].index != -1)
					{
						output_file << problem->x[i][j].index << ":" << problem->x[i][j].value << " " ;					
					}
					++j;
				}
			}
			else
			{
				for(j = 0; j <= problem->l; j++)
				{
					output_file << problem->x[i][j].index << ":" << problem->x[i][j].value << " " ;					
				}
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
		QFileInfo file_info;
		
		file_info.setFile(filename.c_str());
		if (!file_info.exists())
		{
			return NULL;
		}
		if (!file_info.isReadable())
		{
			return NULL;
		}
    if (file_info.size() == 0)
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

	svm_node* LibSVMEncoder::encodeCombinedOligoBordersLibSVMVector(const String& sequence,
																																	const vector<pair<UnsignedInt, UnsignedInt> >& parameters,
																																	const vector<DoubleReal>& sigmas,
																																	const String& allowed_characters,
																																	bool 	strict,
																																	bool length_encoding)
	{
		vector<vector< pair<SignedInt, DoubleReal> > > temp_vectors;
		UnsignedInt                                      number_of_nodes = 0;
		svm_node*                                        nodes = NULL;
		UnsignedInt                                      insertion_index = 0;
		vector<DoubleReal>															 gauss_table;
		svm_node*																				 temp_vector;
				
		temp_vectors.resize(parameters.size(), vector< pair<SignedInt, DoubleReal> >());		
				
		for(UnsignedInt i = 0; i < parameters.size(); ++i)
		{
			if (length_encoding && i == parameters.size() - 1)
			{
				encodeLengthOligo(sequence, temp_vectors[i]);
			}
			else
			{						
				encodeOligoBorders(sequence, parameters[i].first, allowed_characters, parameters[i].second, temp_vectors[i], strict);
			}
			// adding the size of the particular oligo encoding plus one for the end node (-1, 0)
			number_of_nodes += temp_vectors[i].size() + 1;
		}
				
		// reserving the space for the combined vector at the first positions the norm of the particular
		// oligo vectors is stored
		nodes = new svm_node[number_of_nodes + parameters.size()];
		for(insertion_index = 0; insertion_index < parameters.size(); ++insertion_index)
		{
			SVMWrapper::calculateGaussTable(parameters[insertion_index].second, 
																			sigmas[insertion_index], gauss_table);
			temp_vector = encodeLibSVMVector(temp_vectors[insertion_index]);																			
																			
			nodes[insertion_index].index = -1;
			nodes[insertion_index].value = sqrt(SVMWrapper::kernelOligo(temp_vector, temp_vector, gauss_table));
			delete [] temp_vector;
		}
		
		for(UnsignedInt i = 0; i < parameters.size(); ++i)
		{
			for(vector< pair<SignedInt, DoubleReal> >::iterator it = temp_vectors[i].begin();
					it != temp_vectors[i].end();
					++it)
			{
					nodes[insertion_index].index = it->first;
					nodes[insertion_index].value = it->second;
					++insertion_index;
			}
			nodes[insertion_index].index = -1;
			nodes[insertion_index].value = 0;
			++insertion_index;				
		}
		return nodes;
	}

	void LibSVMEncoder::encodeLengthOligo(String sequence,
																	 		  vector< pair<SignedInt, DoubleReal> >& libsvm_vector)
	{
		libsvm_vector.clear();
		libsvm_vector.push_back(make_pair((SignedInt) sequence.size(), 1));		
	}

	void LibSVMEncoder::encodeOligoFeatureVector(const String& 															sequence, 
																							 DoubleReal 																parameter, 
																							 vector<double_pt_2_string_double>  				functions, 
																							 vector< pair<SignedInt, DoubleReal> >& 		encoded_vector, 
																							 UnsignedInt 																start_index,
																							 bool																				length_encoding)
	{
		encoded_vector.clear();
		for(UnsignedInt i = 0; i < functions.size() - 1; ++i)
		{
			encoded_vector.push_back(make_pair(start_index + i, (functions[i](sequence, parameter))));
			encoded_vector.push_back(make_pair(-1, 0));
		}
		encoded_vector.push_back(make_pair(start_index + functions.size() - 1, (functions[functions.size() - 1](sequence, parameter))));
		if (length_encoding && functions.size() > 0)
		{
			encoded_vector.push_back(make_pair(-1, 0));
			encoded_vector.push_back(make_pair(start_index + functions.size(), sequence.size()));
		}
		else if (length_encoding)
		{
			encoded_vector.push_back(make_pair(start_index + functions.size(), sequence.size()));
		}
	}

	void LibSVMEncoder::encodeOligoBorders(String sequence,
													UnsignedInt k_mer_length,
													const String& allowed_characters,
													UnsignedInt border_length,
													vector< pair<SignedInt, DoubleReal> >& libsvm_vector,
													bool strict,
													bool length_encoding)
	{
	  multimap<SignedInt, UnsignedInt>  	         			ordered_tree;
	  multimap<SignedInt, UnsignedInt>::iterator 				elements;
	  pair<SignedInt, UnsignedInt>  	             			values;
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
			  	left_border = (UnsignedInt) (floor(sequence_length - k_mer_length + 1) / 2);
			  	right_border = (UnsignedInt) (ceil(sequence_length - k_mer_length + 1) / 2);
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
		  values.first = -1 * ((SignedInt) (oligo_value + 1));
		  values.second = 1;
		  ordered_tree.insert(values);
		
		  for(UnsignedInt j = 1; j < left_border; j++)
		  {
				oligo_value -= factor * residue_values[sequence[j - 1]];
				oligo_value = oligo_value * number_of_residues + residue_values[sequence[j + k_mer_length - 1]];
		
				values.first = -1 * ((SignedInt) (oligo_value + 1));
				values.second = j + 1;
		
				ordered_tree.insert(values);	
		  }
		  oligo_value = 0;
		  factor = 1;
		  
		  for(UnsignedInt k = right_border + k_mer_length; k > right_border; k--)
		  {
				oligo_value += factor * residue_values[sequence[k - 1]];
				factor *= number_of_residues;
			}
		  factor /= number_of_residues;
		  values.first = oligo_value + 1;
		  values.second = 1;
		  ordered_tree.insert(values);
		  for(UnsignedInt j = right_border + 1; j < sequence_length - k_mer_length + 1; j++)
		  {
				oligo_value -= factor * residue_values[sequence[j - 1]];
				oligo_value = oligo_value * number_of_residues + residue_values[sequence[j + k_mer_length - 1]];
		
				values.first = oligo_value + 1;
				values.second = j - right_border + 1;
		
				ordered_tree.insert(values);	
		  }
		  	
		  for(elements = ordered_tree.begin(); elements != ordered_tree.end(); ++elements)
		  {
				libsvm_vector.push_back(make_pair(elements->second, elements->first));	
		  }
		  if (length_encoding)
		  {
		  	libsvm_vector.push_back(make_pair((SignedInt) sequence.size(), pow(k_mer_length, number_of_residues) + 1));
		  }
		}

	}
	
	svm_problem* LibSVMEncoder::encodeLibSVMProblemWithLengthOligoVectors(const vector<String>&     sequences,
																																				vector<DoubleReal>*       labels)
	{
		vector<svm_node*> vectors;
		vector< pair<SignedInt, DoubleReal> > encoded_vector;
		svm_node* libsvm_vector;
		
		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{			
			encodeLengthOligo(sequences[i], encoded_vector);
			libsvm_vector = encodeLibSVMVector(encoded_vector);
			vectors.push_back(libsvm_vector);
		}
		
		return encodeLibSVMProblem(vectors, labels);				
	}
	
	svm_problem* LibSVMEncoder::encodeLibSVMProblemWithOligoBorderAndFeatureVectors(const vector<String>&     sequences,
																																	 vector<DoubleReal>*		  			labels,
																																	 UnsignedInt 										k_mer_length,
																																	 const String& 	 				  			allowed_characters,
																																	 UnsignedInt 										border_length,
																														 			 vector<double_pt_2_string_double> functions,
																														 			 DoubleReal											ph,
																														 			 bool 													strict,
																														 			 bool 													length_encoding)
	{
		vector<svm_node*> vectors;
		vector< pair<SignedInt, DoubleReal> > encoded_vector1;
		vector< pair<SignedInt, DoubleReal> > encoded_vector2;
		vector< pair<SignedInt, DoubleReal> > encoded_vector3;
		svm_node* libsvm_vector;
		UnsignedInt start_index = 1;
		
		// reserve space for the norms that are stored at the beginning of the libsvm vector
		for(UnsignedInt j = 0; j <= functions.size(); ++j)
		{
			encoded_vector3.push_back(make_pair(-1, 0));
		}

		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{			
			encodeOligoBorders(sequences[i], k_mer_length, allowed_characters, border_length, encoded_vector1, strict);
			if (functions.size() > 0 || length_encoding)
			{
				encoded_vector1.push_back(make_pair(-1, 0));
			}
			encodeOligoFeatureVector(sequences[i], ph, functions, encoded_vector2, start_index, length_encoding);
			encoded_vector1.insert(encoded_vector1.end(), encoded_vector2.begin(), encoded_vector2.end());
			encoded_vector1.insert(encoded_vector1.begin(), encoded_vector3.begin(), encoded_vector3.end());
			libsvm_vector = encodeLibSVMVector(encoded_vector1);
			vectors.push_back(libsvm_vector);
		}
		
		return encodeLibSVMProblem(vectors, labels);				
				
	}																														 			 	

	svm_problem* LibSVMEncoder::encodeLibSVMProblemWithFeatureVectors(const vector<String>&     				sequences,
																																	 vector<DoubleReal>*		  					labels,
																														 			 vector<double_pt_2_string_double> 	functions,
																														 			 DoubleReal													ph,
																														 			 bool 															length_encoding)
	{
		vector<svm_node*> vectors;
		vector< pair<SignedInt, DoubleReal> > encoded_vector1;
		vector< pair<SignedInt, DoubleReal> > encoded_vector2;
		svm_node* libsvm_vector;
		UnsignedInt start_index = 1;
		
		// reserve space for the norms that are stored at the beginning of the libsvm vector
		for(UnsignedInt j = 0; j < functions.size(); ++j)
		{
			encoded_vector2.push_back(make_pair(-1, 0));
		}
		if (length_encoding)
		{
			encoded_vector2.push_back(make_pair(-1, 0));
		}

		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{			
			encodeOligoFeatureVector(sequences[i], ph, functions, encoded_vector1, start_index, length_encoding);
			encoded_vector1.insert(encoded_vector1.begin(), encoded_vector2.begin(), encoded_vector2.end());
			libsvm_vector = encodeLibSVMVector(encoded_vector1);
						
			vectors.push_back(libsvm_vector);
		}
		
		return encodeLibSVMProblem(vectors, labels);				
				
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
		
		while (vector[i].index != -1)
		{
			output = output + "(" + String(vector[i].index) + ", " + String(vector[i].value) + ") ";
			++i;
		}
	}
	
	void LibSVMEncoder::oligoBorderVectorToString(svm_node* vector, UnsignedInt border_length, String& output)
	{
		map<SignedInt, DoubleReal> left_part;
		map<SignedInt, DoubleReal> right_part;
		
		UnsignedInt i = 0;
		UnsignedInt zero_counter = 0;
		
		output = "";
		if (vector != NULL)
		{

			while(vector[i].index != -1 && vector[i].value < 0)
			{
				left_part.insert(make_pair(vector[i].index, -1 * vector[i].value));
				++i;
			}

			while(vector[i].index != -1)
			{
				right_part.insert(make_pair(vector[i].index, vector[i].value));
				++i;
			}
			i = 0;
			for(map<SignedInt, DoubleReal>::iterator it = left_part.begin();
					it != left_part.end();
					++it)
			{
				output = output + String(it->second) + " ";
				++i;
			}
			while(i < border_length)
			{
				output = output + "0 ";
				++zero_counter;
				++i;
			}
			while(zero_counter > 0)
			{
				output = output + "0 ";
				--zero_counter;
			}
			for(map<SignedInt, DoubleReal>::iterator it = right_part.begin();
					it != right_part.end();
					++it)
			{
				output = output + String(it->second) + " ";
			}
		}
	}
	
	void LibSVMEncoder::combinedOligoBorderVectorToString(svm_node* vector, UnsignedInt number_of_combinations, String& output)
	{
		map<SignedInt, DoubleReal> left_part;
		map<SignedInt, DoubleReal> right_part;
		
		UnsignedInt i = 0;
		UnsignedInt end_counter = 0;
		
		output = "";
	
		if (vector != NULL)
		{
			while(i < number_of_combinations)
			{
				output = output + "(" + String(vector[i].index) + "," + String(vector[i].value) + ") ";				
				++i;
			}

			while(end_counter < number_of_combinations)
			{
				if (vector[i].index == -1)
				{
					++end_counter;
				}
				output = output + "(" + String(vector[i].index) + "," + String(vector[i].value) + ") ";
				++i;
			}
		}
	}
	
	void LibSVMEncoder::libSVMVectorsToString(svm_problem* vector, String& output)
	{
		String temp_string = "";
		
		output = "";
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
	
	void LibSVMEncoder::oligoBorderVectorsToString(svm_problem* vector, UnsignedInt border_length, String& output)
	{
		String temp_string = "";
		
		output = "";
		if (vector != NULL)
		{
			for(SignedInt i = 0; i < vector->l; ++i)
			{
				oligoBorderVectorToString(vector->x[i], border_length, temp_string);
				output = output + temp_string + "\n";
			}
		}
	}
	
	void LibSVMEncoder::combinedOligoBorderVectorsToString(svm_problem* vector, UnsignedInt number_of_combinations, String& output)
	{
		String temp_string = "";
		
		output = "";
		if (vector != NULL)
		{
			for(SignedInt i = 0; i < vector->l; ++i)
			{
				combinedOligoBorderVectorToString(vector->x[i], number_of_combinations, temp_string);
				output = output + temp_string + "\n";
			}
		}
	}
	
	svm_problem* LibSVMEncoder::encodeLibSVMProblemWithCombinedOligoBorderVectors(const vector<String>&     sequences,
																																 vector<DoubleReal>*  		 labels,
																																 const vector<pair<UnsignedInt, UnsignedInt> >& parameters,
																																 const vector<DoubleReal>& sigmas,
																																 const String& 	 				   allowed_characters,
																																 bool											 strict,
																																 bool 										 length_encoding)	
	{
		vector<svm_node*> vectors;
		svm_node* libsvm_vector;
		
		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{			
			libsvm_vector = encodeCombinedOligoBordersLibSVMVector(sequences[i], parameters, sigmas, allowed_characters, strict, length_encoding);
			vectors.push_back(libsvm_vector);
		}
		
		return encodeLibSVMProblem(vectors, labels);		
		
	}

	svm_problem* LibSVMEncoder::encodeLibSVMProblemWithCompositionLengthAndHydroVectors(const vector<String>&    sequences,
																																								 std::vector<DoubleReal>* labels,
																																								 const String& 	 				  allowed_characters,
																																								 UnsignedInt 							maximum_sequence_length)
	{
		vector<svm_node*> vectors;
		vector< pair<SignedInt, DoubleReal> >* encoded_vector;
		svm_node* libsvm_vector;
		DoubleReal 							sum1 								= 0;
		DoubleReal 							sum2 								= 0;
		DoubleReal 							pi	 								= 3.1415926536;
		map<char, DoubleReal> 	hydrophobicities;
		
		hydrophobicities.insert(make_pair('A', 0.61));
		hydrophobicities.insert(make_pair('L', 1.53));
		hydrophobicities.insert(make_pair('R', 0.60));
		hydrophobicities.insert(make_pair('K', 1.15));
		hydrophobicities.insert(make_pair('N', 0.06));
		hydrophobicities.insert(make_pair('M', 1.18));
		hydrophobicities.insert(make_pair('D', 0.46));
		hydrophobicities.insert(make_pair('F', 2.02));
		hydrophobicities.insert(make_pair('C', 1.07));
		hydrophobicities.insert(make_pair('P', 1.95));
		hydrophobicities.insert(make_pair('Q', 0.));
		hydrophobicities.insert(make_pair('S', 0.05));
		hydrophobicities.insert(make_pair('E', 0.47));
		hydrophobicities.insert(make_pair('T', 0.05));
		hydrophobicities.insert(make_pair('G', 0.07));
		hydrophobicities.insert(make_pair('W', 2.65));
		hydrophobicities.insert(make_pair('H', 0.61));
		hydrophobicities.insert(make_pair('Y', 1.88));
		hydrophobicities.insert(make_pair('I', 2.22));
		hydrophobicities.insert(make_pair('V', 1.32));

		for(UnsignedInt i = 0; i < sequences.size(); i++)
		{
			sum1 = 0;
			sum2 = 0;
			for(UnsignedInt j = 0; j < sequences[i].size(); ++j)
			{
				sum1 += hydrophobicities[sequences[i].at(j)] * sin(2 * (j + 1) * pi / 3.6);
				sum2 += hydrophobicities[sequences[i].at(j)] * cos(2 * (j + 1) * pi / 3.6);
			}
			sum1 *= sum1;
			sum2 *= sum2;
			
			encoded_vector = encodeCompositionVector(sequences[i], allowed_characters);
			encoded_vector->push_back(make_pair(allowed_characters.size() + 1, ((DoubleReal) sequences[i].length()) / maximum_sequence_length));
			encoded_vector->push_back(make_pair(allowed_characters.size() + 2, sqrt(sum1 + sum2)));
			libsvm_vector = encodeLibSVMVector(*encoded_vector);
			vectors.push_back(libsvm_vector);
		}
		
		return encodeLibSVMProblem(vectors, labels);		
	}
																																		 	
} // namespace OpenMS
