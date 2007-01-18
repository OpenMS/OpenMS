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

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <svm.h>

///////////////////////////

#include <string>
#include <vector>

///////////////////////////

START_TEST(LibSVMEncoder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

LibSVMEncoder* ptr;
LibSVMEncoder encoder;

CHECK((LibSVMEncoder()))
  ptr = new LibSVMEncoder();
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~LibSVMEncoder()))
	delete ptr;
RESULT

CHECK((std::vector< std::pair<SignedInt, DoubleReal>* encodeCompositionVector(const String& sequence, const String& allowed_characters = "ACDEFGHIKLMNPQRSTVWY")))
	String sequence = "ACCGGGTTTT";
	String allowed_characters = "ACNGT";
	vector< pair<SignedInt, DoubleReal> >* encoded_sequence;
	std::vector< std::pair<SignedInt, DoubleReal> >::iterator it;
			
	encoded_sequence = encoder.encodeCompositionVector(sequence, allowed_characters);
	it = encoded_sequence->begin(); 
	TEST_EQUAL(it->first, 1)
	TEST_REAL_EQUAL(it->second, 0.1)
	it++;
	TEST_EQUAL(it == encoded_sequence->end(), false)
	TEST_EQUAL(it->first, 2)
	TEST_REAL_EQUAL(it->second, 0.2)
	it++;
	TEST_EQUAL(it == encoded_sequence->end(), false)
	TEST_EQUAL(it->first, 4)
	TEST_REAL_EQUAL(it->second, 0.3)
	it++;
	TEST_EQUAL(it == encoded_sequence->end(), false)
	TEST_EQUAL(it->first, 5)
	TEST_REAL_EQUAL(it->second, 0.4)
	it++;
	TEST_EQUAL(it == encoded_sequence->end(), true)
	delete encoded_sequence;
RESULT

CHECK((std::vector< std::vector< std::pair<SignedInt, DoubleReal>* encodeCompositionVectors(const std::vector<String>& sequences, const String& allowed_characters)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector< vector< pair<SignedInt, DoubleReal> > >* encoded_sequences;
	vector< pair<SignedInt, DoubleReal> >::iterator it;
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoded_sequences = encoder.encodeCompositionVectors(sequences, allowed_characters);
	it = (*encoded_sequences)[0].begin(); 
	TEST_EQUAL(it->first, 1)
	TEST_REAL_EQUAL(it->second, 0.1)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[0].end(), false)
	TEST_EQUAL(it->first, 2)
	TEST_REAL_EQUAL(it->second, 0.2)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[0].end(), false)
	TEST_EQUAL(it->first, 4)
	TEST_REAL_EQUAL(it->second, 0.3)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[0].end(), false)
	TEST_EQUAL(it->first, 5)
	TEST_REAL_EQUAL(it->second, 0.4)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[0].end(), true)
	it = (*encoded_sequences)[1].begin(); 
	TEST_EQUAL(it == (*encoded_sequences)[1].end(), false)
	TEST_EQUAL(it->first, 1)
	TEST_REAL_EQUAL(it->second, 0.5)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[1].end(), false)
	TEST_EQUAL(it->first, 2)
	TEST_REAL_EQUAL(it->second, 0.5)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[1].end(), true)
	delete encoded_sequences;
RESULT

CHECK((std::vector<svm_node*>* encodeLibSVMVectors( const std::vector< std::vector< std::pair<SignedInt, DoubleReal> > >& feature_vectors)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector<vector< pair<SignedInt, DoubleReal> > >* encoded_sequences;
	vector<svm_node*>* libsvm_sequences;
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoded_sequences = encoder.encodeCompositionVectors(sequences, allowed_characters);
	libsvm_sequences = encoder.encodeLibSVMVectors(*encoded_sequences);
	nodes = (*libsvm_sequences)[0];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_EQUAL(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_EQUAL(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_EQUAL(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_EQUAL(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	nodes = (*libsvm_sequences)[1];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_EQUAL(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_EQUAL(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, -1)
	delete encoded_sequences;
	delete libsvm_sequences;
	delete nodes;
	
RESULT

CHECK((svm_node* encodeLibSVMVector( const std::vector< std::pair<SignedInt, DoubleReal> >& feature_vector)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector< pair<SignedInt, DoubleReal> >* encoded_sequence;
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoded_sequence = encoder.encodeCompositionVector(sequences[0], allowed_characters);
	nodes = encoder.encodeLibSVMVector(*encoded_sequence);
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_EQUAL(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_EQUAL(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_EQUAL(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_EQUAL(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	
	delete encoded_sequence;
	delete nodes;
	
RESULT

CHECK((svm_problem* encodeLibSVMProblem(const std::vector<svm_node*>& vectors, std::vector<DoubleReal>* labels)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector<vector< pair<SignedInt, DoubleReal> > >* encoded_sequences;
	vector<svm_node*>* libsvm_sequences;
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	svm_problem* problem;
	vector<DoubleReal>* labels = new vector<DoubleReal>();
	
	labels->push_back(2.1);
	labels->push_back(1.3);

	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoded_sequences = encoder.encodeCompositionVectors(sequences, allowed_characters);
	libsvm_sequences = encoder.encodeLibSVMVectors(*encoded_sequences);
	problem = encoder.encodeLibSVMProblem(*libsvm_sequences, labels);
	TEST_EQUAL(problem->l, 2)
	TEST_REAL_EQUAL(problem->y[0], 2.1);
	TEST_REAL_EQUAL(problem->y[1], 1.3);
	nodes = problem->x[0];
	TEST_REAL_EQUAL(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_EQUAL(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_EQUAL(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_EQUAL(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	nodes = problem->x[1];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_EQUAL(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_EQUAL(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, -1)
	delete encoded_sequences;
	delete libsvm_sequences;
	delete nodes;
	delete problem;

RESULT

CHECK((svm_problem* encodeLibSVMProblemWithCompositionAndLengthVectors(const std::vector<String>& sequences, std::vector<DoubleReal>* labels, const String& allowed_characters, UnsignedInt maximum_sequence_length)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	svm_problem* problem;
	vector<DoubleReal>* labels = new vector<DoubleReal>();
	
	labels->push_back(2.1);
	labels->push_back(1.3);

	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	problem = encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(sequences, labels, allowed_characters, 10);
	TEST_EQUAL(problem->l, 2)
	TEST_REAL_EQUAL(problem->y[0], 2.1);
	TEST_REAL_EQUAL(problem->y[1], 1.3);
	nodes = problem->x[0];
	TEST_REAL_EQUAL(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_EQUAL(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_EQUAL(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_EQUAL(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, 6)
	TEST_REAL_EQUAL(nodes[4].value, 1)
	TEST_EQUAL(nodes[5].index, -1)
	nodes = problem->x[1];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_EQUAL(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_EQUAL(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, 6)
	TEST_REAL_EQUAL(nodes[2].value, 0.4)
	TEST_EQUAL(nodes[3].index, -1)
	delete labels;
	delete nodes;
	delete problem;
RESULT

CHECK((svm_problem* encodeLibSVMProblemWithCompositionVectors(const std::vector<String>& sequences, std::vector<DoubleReal>* labels, const String& allowed_characters)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	svm_problem* problem;
	vector<DoubleReal>* labels = new vector<DoubleReal>();
	
	labels->push_back(2.1);
	labels->push_back(1.3);

	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	problem = encoder.encodeLibSVMProblemWithCompositionVectors(sequences, labels, allowed_characters);
	TEST_EQUAL(problem->l, 2)
	TEST_REAL_EQUAL(problem->y[0], 2.1);
	TEST_REAL_EQUAL(problem->y[1], 1.3);
	nodes = problem->x[0];
	TEST_REAL_EQUAL(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_EQUAL(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_EQUAL(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_EQUAL(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	nodes = problem->x[1];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_EQUAL(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_EQUAL(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, -1)
	delete labels;
	delete problem;
RESULT

CHECK((bool storeLibSVMProblem(const String& filename, const svm_problem* problem) const))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	svm_problem* problem;
	vector<DoubleReal>* labels = new vector<DoubleReal>();
	String temp_filename = "data/LibSVMEncoder_test.tmp";
	
	labels->push_back(2.1);
	labels->push_back(1.3);

	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	problem = encoder.encodeLibSVMProblemWithCompositionVectors(sequences, labels, allowed_characters);
	NEW_TMP_FILE(temp_filename)
	encoder.storeLibSVMProblem(temp_filename, problem);
	TEST_FILE("data/LibSVMEncoder_test.txt", temp_filename.c_str())

RESULT

CHECK(svm_problem* loadLibSVMProblem(const String& filename))
	String allowed_characters = "ACNGT";
	svm_problem* problem;
	String temp_filename = "data/LibSVMEncoder_test.tmp";
	
	NEW_TMP_FILE(temp_filename)
	problem = encoder.loadLibSVMProblem("data/LibSVMEncoder_test.txt");
	encoder.storeLibSVMProblem(temp_filename, problem);
	TEST_FILE("data/LibSVMEncoder_test.txt", temp_filename.c_str())

RESULT

CHECK((void encodeOligoBorders(String sequence, UnsignedInt k_mer_length, const String& allowed_characters, UnsignedInt border_length, std::vector< std::pair<SignedInt, DoubleReal> >& libsvm_vector, bool strict = false, bool length_encoding = false)))
	String sequence = "ACNNGTATCA";
	String allowed_characters = "ACNGT";
	String output;
	UnsignedInt border_length = 3;
	vector< pair<SignedInt, DoubleReal> > encoded_sequence;
	
	encoder.encodeOligoBorders(sequence, 1, allowed_characters, border_length, encoded_sequence);
	encoder.libSVMVectorToString(encoder.encodeLibSVMVector(encoded_sequence), output);
	TEST_EQUAL(output, "(2, 1) (2, 1) (3, 2) (3, 2) (4, 3) (6, 3) ")
	encoder.encodeOligoBorders(sequence, 2, allowed_characters, border_length, encoded_sequence);
	encoder.libSVMVectorToString(encoder.encodeLibSVMVector(encoded_sequence), output);
	TEST_EQUAL(output, "(3, 1) (3, 1) (9, 2) (11, 2) (14, 3) (22, 3) ")
RESULT

CHECK((svm_problem* encodeLibSVMProblemWithOligoBorderVectors(const std::vector<String>& sequences, std::vector<DoubleReal>* labels, UnsignedInt k_mer_length, const String& allowed_characters, UnsignedInt border_length, bool strict = false, bool length_encoding = false)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	String output;
	UnsignedInt border_length = 3;
	vector< pair<SignedInt, DoubleReal> > encoded_sequence;
  vector<DoubleReal> labels;
  struct svm_problem* data;

  labels.push_back(1);
  labels.push_back(2);
  sequences.push_back("ACNNGTATCA");
  sequences.push_back("AACNNGTACCA");
	data = encoder.encodeLibSVMProblemWithOligoBorderVectors(sequences, &labels, 1, allowed_characters, border_length);
	encoder.libSVMVectorToString(data->x[0], output);
	TEST_EQUAL(output, "(2, 1) (2, 1) (3, 2) (3, 2) (4, 3) (6, 3) ")
	encoder.libSVMVectorToString(data->x[1], output);
	TEST_EQUAL(output, "(2, 1) (2, 2) (2, 1) (3, 3) (3, 3) (3, 2) ")
RESULT

CHECK((void libSVMVectorToString(svm_node* vector, String& output)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector< pair<SignedInt, DoubleReal> >* encoded_sequence;
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	String output;	
	String correct_output = "(1, 0.1) (2, 0.2) (4, 0.3) (5, 0.4) ";
	
	sequences.push_back(String("ACCGGGTTTT"));			
			
	encoded_sequence = encoder.encodeCompositionVector(sequences[0], allowed_characters);
	nodes = encoder.encodeLibSVMVector(*encoded_sequence);
	
	encoder.libSVMVectorToString(nodes, output);
	
	TEST_EQUAL(output, correct_output)
	
RESULT

CHECK((void libSVMVectorsToString(svm_problem* vector, String& output)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	String output;	
	String correct_output = "(1, 0.1) (2, 0.2) (4, 0.3) (5, 0.4) \n(1, 0.5) (2, 0.5) \n";	
	vector<DoubleReal>* labels = new vector<DoubleReal>();
	svm_problem* problem;
	
	labels->push_back(2.1);
	labels->push_back(1.3);
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));
		
	problem = encoder.encodeLibSVMProblemWithCompositionVectors(sequences, labels, allowed_characters);			
	encoder.libSVMVectorsToString(problem, output);
	TEST_EQUAL(output, correct_output)	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
