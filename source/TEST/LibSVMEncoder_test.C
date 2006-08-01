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

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <svm.h>

///////////////////////////

#include <string>
#include <vector>

///////////////////////////

START_TEST(LibSVMEncoder, "$Id: LibSVMEncoder_test.C 200 2006-07-24 14:14:07Z nicopfeifer $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

LibSVMEncoder* ptr;
LibSVMEncoder encoder;

CHECK(LibSVMEncoder())
  ptr = new LibSVMEncoder();
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~LibSVMEncoder())
	delete ptr;
RESULT

CHECK((std::vector< std::pair<UnsignedInt, DoubleReal>* encodeCompositionVector(const String& sequence, const String& allowed_characters)))
	String sequence = "ACCGGGTTTT";
	String allowed_characters = "ACNGT";
	vector< pair<UnsignedInt, DoubleReal> >* encoded_sequence;
	std::vector< std::pair<UnsignedInt, DoubleReal> >::iterator it;
			
	encoded_sequence = encoder.encodeCompositionVector(sequence, allowed_characters);
	it = encoded_sequence->begin(); 
	TEST_EQUAL(it->first, 0)
	TEST_REAL_EQUAL(it->second, 0.1)
	it++;
	TEST_EQUAL(it == encoded_sequence->end(), false)
	TEST_EQUAL(it->first, 1)
	TEST_REAL_EQUAL(it->second, 0.2)
	it++;
	TEST_EQUAL(it == encoded_sequence->end(), false)
	TEST_EQUAL(it->first, 3)
	TEST_REAL_EQUAL(it->second, 0.3)
	it++;
	TEST_EQUAL(it == encoded_sequence->end(), false)
	TEST_EQUAL(it->first, 4)
	TEST_REAL_EQUAL(it->second, 0.4)
	it++;
	TEST_EQUAL(it == encoded_sequence->end(), true)
	delete encoded_sequence;
RESULT

CHECK((std::vector< std::vector< std::pair<UnsignedInt, DoubleReal>* encodeCompositionVectors(const std::vector<String>& sequences, const String& allowed_characters)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector< vector< pair<UnsignedInt, DoubleReal> > >* encoded_sequences;
	vector< pair<UnsignedInt, DoubleReal> >::iterator it;
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoded_sequences = encoder.encodeCompositionVectors(sequences, allowed_characters);
	it = (*encoded_sequences)[0].begin(); 
	TEST_EQUAL(it->first, 0)
	TEST_REAL_EQUAL(it->second, 0.1)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[0].end(), false)
	TEST_EQUAL(it->first, 1)
	TEST_REAL_EQUAL(it->second, 0.2)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[0].end(), false)
	TEST_EQUAL(it->first, 3)
	TEST_REAL_EQUAL(it->second, 0.3)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[0].end(), false)
	TEST_EQUAL(it->first, 4)
	TEST_REAL_EQUAL(it->second, 0.4)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[0].end(), true)
	it = (*encoded_sequences)[1].begin(); 
	TEST_EQUAL(it == (*encoded_sequences)[1].end(), false)
	TEST_EQUAL(it->first, 0)
	TEST_REAL_EQUAL(it->second, 0.5)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[1].end(), false)
	TEST_EQUAL(it->first, 1)
	TEST_REAL_EQUAL(it->second, 0.5)
	it++;
	TEST_EQUAL(it == (*encoded_sequences)[1].end(), true)
	delete encoded_sequences;
RESULT

CHECK((std::vector<svm_node*>* encodeLIBSVMVectors( const std::vector< std::vector< std::pair<UnsignedInt, DoubleReal> > >& feature_vectors)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector<vector< pair<UnsignedInt, DoubleReal> > >* encoded_sequences;
	vector<svm_node*>* libsvm_sequences;
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoded_sequences = encoder.encodeCompositionVectors(sequences, allowed_characters);
	libsvm_sequences = encoder.encodeLIBSVMVectors(*encoded_sequences);
	nodes = (*libsvm_sequences)[0];
	TEST_EQUAL(nodes[0].index, 0)
	TEST_REAL_EQUAL(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 1)
	TEST_REAL_EQUAL(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 3)
	TEST_REAL_EQUAL(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 4)
	TEST_REAL_EQUAL(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	nodes = (*libsvm_sequences)[1];
	TEST_EQUAL(nodes[0].index, 0)
	TEST_REAL_EQUAL(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 1)
	TEST_REAL_EQUAL(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, -1)
	delete encoded_sequences;
	delete libsvm_sequences;
	delete nodes;
	
RESULT

CHECK((svm_node* encodeLIBSVMVector( const std::vector< std::pair<UnsignedInt, DoubleReal> >& feature_vector)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector< pair<UnsignedInt, DoubleReal> >* encoded_sequence;
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoded_sequence = encoder.encodeCompositionVector(sequences[0], allowed_characters);
	nodes = encoder.encodeLIBSVMVector(*encoded_sequence);
	TEST_EQUAL(nodes[0].index, 0)
	TEST_REAL_EQUAL(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 1)
	TEST_REAL_EQUAL(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 3)
	TEST_REAL_EQUAL(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 4)
	TEST_REAL_EQUAL(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	
	delete encoded_sequence;
	delete nodes;
	
RESULT

CHECK((svm_problem* encodeLIBSVMProblem(const std::vector<svm_node*>& vectors, std::vector<DoubleReal>* labels)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector<vector< pair<UnsignedInt, DoubleReal> > >* encoded_sequences;
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
	libsvm_sequences = encoder.encodeLIBSVMVectors(*encoded_sequences);
	problem = encoder.encodeLIBSVMProblem(*libsvm_sequences, labels);
	TEST_EQUAL(problem->l, 2)
	TEST_REAL_EQUAL(problem->y[0], 2.1);
	TEST_REAL_EQUAL(problem->y[1], 1.3);
	nodes = problem->x[0];
	TEST_REAL_EQUAL(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 1)
	TEST_REAL_EQUAL(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 3)
	TEST_REAL_EQUAL(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 4)
	TEST_REAL_EQUAL(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	nodes = problem->x[1];
	TEST_EQUAL(nodes[0].index, 0)
	TEST_REAL_EQUAL(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 1)
	TEST_REAL_EQUAL(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, -1)
	delete encoded_sequences;
	delete libsvm_sequences;
	delete nodes;
	delete problem;

RESULT

CHECK((svm_problem* encodeLIBSVMProblemWithCompositionAndLengthVectors(const std::vector<String>& sequences, std::vector<DoubleReal>* labels, const String& allowed_characters, UnsignedInt maximum_sequence_length)))
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
			
	problem = encoder.encodeLIBSVMProblemWithCompositionAndLengthVectors(sequences, labels, allowed_characters, 10);
	TEST_EQUAL(problem->l, 2)
	TEST_REAL_EQUAL(problem->y[0], 2.1);
	TEST_REAL_EQUAL(problem->y[1], 1.3);
	nodes = problem->x[0];
	TEST_REAL_EQUAL(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 1)
	TEST_REAL_EQUAL(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 3)
	TEST_REAL_EQUAL(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 4)
	TEST_REAL_EQUAL(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, 5)
	TEST_REAL_EQUAL(nodes[4].value, 1)
	TEST_EQUAL(nodes[5].index, -1)
	nodes = problem->x[1];
	TEST_EQUAL(nodes[0].index, 0)
	TEST_REAL_EQUAL(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 1)
	TEST_REAL_EQUAL(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, 5)
	TEST_REAL_EQUAL(nodes[2].value, 0.4)
	TEST_EQUAL(nodes[3].index, -1)
	delete labels;
	delete nodes;
	delete problem;
RESULT

CHECK((svm_problem* encodeLIBSVMProblemWithCompositionVectors(const std::vector<String>& sequences, std::vector<DoubleReal>* labels, const String& allowed_characters)))
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
			
	problem = encoder.encodeLIBSVMProblemWithCompositionVectors(sequences, labels, allowed_characters);
	TEST_EQUAL(problem->l, 2)
	TEST_REAL_EQUAL(problem->y[0], 2.1);
	TEST_REAL_EQUAL(problem->y[1], 1.3);
	nodes = problem->x[0];
	TEST_REAL_EQUAL(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 1)
	TEST_REAL_EQUAL(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 3)
	TEST_REAL_EQUAL(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 4)
	TEST_REAL_EQUAL(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	nodes = problem->x[1];
	TEST_EQUAL(nodes[0].index, 0)
	TEST_REAL_EQUAL(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 1)
	TEST_REAL_EQUAL(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, -1)
	delete labels;
	delete problem;
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
