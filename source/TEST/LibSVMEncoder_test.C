// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

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

START_SECTION((LibSVMEncoder()))
  ptr = new LibSVMEncoder();
  TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~LibSVMEncoder()))
	delete ptr;
END_SECTION

START_SECTION((void encodeCompositionVector(const String &sequence, std::vector< std::pair< Int, DoubleReal > > &encoded_vector, const String &allowed_characters="ACDEFGHIKLMNPQRSTVWY")))
	String sequence = "ACCGGGTTTT";
	String allowed_characters = "ACNGT";
	vector< pair<Int, DoubleReal> > encoded_sequence;
	std::vector< std::pair<Int, DoubleReal> >::iterator it;
			
	encoder.encodeCompositionVector(sequence, encoded_sequence, allowed_characters);
	it = encoded_sequence.begin(); 
	TEST_EQUAL(it->first, 1)
	TEST_REAL_SIMILAR(it->second, 0.1)
	it++;
	TEST_EQUAL(it == encoded_sequence.end(), false)
	TEST_EQUAL(it->first, 2)
	TEST_REAL_SIMILAR(it->second, 0.2)
	it++;
	TEST_EQUAL(it == encoded_sequence.end(), false)
	TEST_EQUAL(it->first, 4)
	TEST_REAL_SIMILAR(it->second, 0.3)
	it++;
	TEST_EQUAL(it == encoded_sequence.end(), false)
	TEST_EQUAL(it->first, 5)
	TEST_REAL_SIMILAR(it->second, 0.4)
	it++;
	TEST_EQUAL(it == encoded_sequence.end(), true)
END_SECTION

START_SECTION((void encodeCompositionVectors(const std::vector< String > &sequences, const String &allowed_characters, std::vector< std::vector< std::pair< Int, DoubleReal > > > &composition_vectors)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector< vector< pair<Int, DoubleReal> > > encoded_sequences;
	vector< pair<Int, DoubleReal> >::iterator it;
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoder.encodeCompositionVectors(sequences, allowed_characters, encoded_sequences);
	it = encoded_sequences[0].begin(); 
	TEST_EQUAL(it->first, 1)
	TEST_REAL_SIMILAR(it->second, 0.1)
	it++;
	TEST_EQUAL(it == encoded_sequences[0].end(), false)
	TEST_EQUAL(it->first, 2)
	TEST_REAL_SIMILAR(it->second, 0.2)
	it++;
	TEST_EQUAL(it == encoded_sequences[0].end(), false)
	TEST_EQUAL(it->first, 4)
	TEST_REAL_SIMILAR(it->second, 0.3)
	it++;
	TEST_EQUAL(it == encoded_sequences[0].end(), false)
	TEST_EQUAL(it->first, 5)
	TEST_REAL_SIMILAR(it->second, 0.4)
	it++;
	TEST_EQUAL(it == encoded_sequences[0].end(), true)
	it = encoded_sequences[1].begin(); 
	TEST_EQUAL(it == encoded_sequences[1].end(), false)
	TEST_EQUAL(it->first, 1)
	TEST_REAL_SIMILAR(it->second, 0.5)
	it++;
	TEST_EQUAL(it == encoded_sequences[1].end(), false)
	TEST_EQUAL(it->first, 2)
	TEST_REAL_SIMILAR(it->second, 0.5)
	it++;
	TEST_EQUAL(it == encoded_sequences[1].end(), true)
END_SECTION

START_SECTION((void encodeLibSVMVectors(const std::vector< std::vector< std::pair< Int, DoubleReal > > > &feature_vectors, std::vector< svm_node * > &libsvm_vectors)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector<vector< pair<Int, DoubleReal> > > encoded_sequences;
	vector<svm_node*> libsvm_sequences;
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoder.encodeCompositionVectors(sequences, allowed_characters, encoded_sequences);
	encoder.encodeLibSVMVectors(encoded_sequences, libsvm_sequences);
	nodes = libsvm_sequences[0];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_SIMILAR(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_SIMILAR(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_SIMILAR(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	nodes = libsvm_sequences[1];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_SIMILAR(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, -1)
	delete nodes;
	
END_SECTION

START_SECTION((svm_node* encodeLibSVMVector( const std::vector< std::pair<Int, DoubleReal> >& feature_vector)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector< pair<Int, DoubleReal> > encoded_sequence;
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoder.encodeCompositionVector(sequences[0], encoded_sequence, allowed_characters);
	nodes = encoder.encodeLibSVMVector(encoded_sequence);
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_SIMILAR(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_SIMILAR(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_SIMILAR(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	
	delete nodes;
	
END_SECTION

START_SECTION((svm_problem* encodeLibSVMProblem(const std::vector< svm_node * > &vectors, std::vector< DoubleReal > &labels)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector<vector< pair<Int, DoubleReal> > > encoded_sequences;
	vector<svm_node*> libsvm_sequences;
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	svm_problem* problem;
	vector<DoubleReal> labels;
	
	labels.push_back(2.1);
	labels.push_back(1.3);

	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	encoder.encodeCompositionVectors(sequences, allowed_characters, encoded_sequences);
	encoder.encodeLibSVMVectors(encoded_sequences, libsvm_sequences);
	problem = encoder.encodeLibSVMProblem(libsvm_sequences, labels);
	TEST_EQUAL(problem->l, 2)
	TEST_REAL_SIMILAR(problem->y[0], 2.1);
	TEST_REAL_SIMILAR(problem->y[1], 1.3);
	nodes = problem->x[0];
	TEST_REAL_SIMILAR(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_SIMILAR(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_SIMILAR(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	nodes = problem->x[1];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_SIMILAR(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, -1)
	delete nodes;
	delete problem;

END_SECTION

START_SECTION((svm_problem* encodeLibSVMProblemWithCompositionAndLengthVectors(const std::vector< String > &sequences, std::vector< DoubleReal > &labels, const String &allowed_characters, UInt maximum_sequence_length)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	svm_problem* problem;
	vector<DoubleReal> labels;
	
	labels.push_back(2.1);
	labels.push_back(1.3);

	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	problem = encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(sequences, labels, allowed_characters, 10);
	TEST_EQUAL(problem->l, 2)
	TEST_REAL_SIMILAR(problem->y[0], 2.1);
	TEST_REAL_SIMILAR(problem->y[1], 1.3);
	nodes = problem->x[0];
	TEST_REAL_SIMILAR(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_SIMILAR(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_SIMILAR(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, 6)
	TEST_REAL_SIMILAR(nodes[4].value, 1)
	TEST_EQUAL(nodes[5].index, -1)
	nodes = problem->x[1];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_SIMILAR(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, 6)
	TEST_REAL_SIMILAR(nodes[2].value, 0.4)
	TEST_EQUAL(nodes[3].index, -1)
	delete nodes;
	delete problem;
END_SECTION

START_SECTION((svm_problem* encodeLibSVMProblemWithCompositionLengthAndWeightVectors(const std::vector< String > &sequences, std::vector< DoubleReal > &labels, const String &allowed_characters)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	svm_problem* problem;
	vector<DoubleReal> labels;
	
	labels.push_back(2.1);
	labels.push_back(1.3);

	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	problem = encoder.encodeLibSVMProblemWithCompositionLengthAndWeightVectors(sequences, labels, allowed_characters);
	TEST_EQUAL(problem->l, 2)
	TEST_REAL_SIMILAR(problem->y[0], 2.1);
	TEST_REAL_SIMILAR(problem->y[1], 1.3);
	nodes = problem->x[0];
	TEST_REAL_SIMILAR(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_SIMILAR(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_SIMILAR(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, 6)
	TEST_REAL_SIMILAR(nodes[4].value, 10)
	TEST_EQUAL(nodes[5].index, 7)
	TEST_REAL_SIMILAR(nodes[5].value, 870.948464870453)
	TEST_EQUAL(nodes[6].index, -1)
	nodes = problem->x[1];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_SIMILAR(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, 6)
	TEST_REAL_SIMILAR(nodes[2].value, 4)
	TEST_EQUAL(nodes[3].index, 7)
	TEST_REAL_SIMILAR(nodes[3].value, 366.45688)
	TEST_EQUAL(nodes[4].index, -1)
	delete nodes;
	delete problem;
END_SECTION

START_SECTION((svm_problem* encodeLibSVMProblemWithCompositionVectors(const std::vector< String > &sequences, std::vector< DoubleReal > &labels, const String &allowed_characters)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	svm_problem* problem;
	vector<DoubleReal> labels;
	
	labels.push_back(2.1);
	labels.push_back(1.3);

	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	problem = encoder.encodeLibSVMProblemWithCompositionVectors(sequences, labels, allowed_characters);
	TEST_EQUAL(problem->l, 2)
	TEST_REAL_SIMILAR(problem->y[0], 2.1);
	TEST_REAL_SIMILAR(problem->y[1], 1.3);
	nodes = problem->x[0];
	TEST_REAL_SIMILAR(nodes[0].value, 0.1)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.2)
	TEST_EQUAL(nodes[2].index, 4)
	TEST_REAL_SIMILAR(nodes[2].value, 0.3)
	TEST_EQUAL(nodes[3].index, 5)
	TEST_REAL_SIMILAR(nodes[3].value, 0.4)
	TEST_EQUAL(nodes[4].index, -1)
	nodes = problem->x[1];
	TEST_EQUAL(nodes[0].index, 1)
	TEST_REAL_SIMILAR(nodes[0].value, 0.5)
	TEST_EQUAL(nodes[1].index, 2)
	TEST_REAL_SIMILAR(nodes[1].value, 0.5)
	TEST_EQUAL(nodes[2].index, -1)
	delete problem;
END_SECTION

START_SECTION((bool storeLibSVMProblem(const String& filename, const svm_problem* problem) const))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	svm_problem* problem;
	vector<DoubleReal> labels;
	String temp_filename = OPENMS_GET_TEST_DATA_PATH("LibSVMEncoder_test.tmp");
	
	labels.push_back(2.1);
	labels.push_back(1.3);

	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));			
			
	problem = encoder.encodeLibSVMProblemWithCompositionVectors(sequences, labels, allowed_characters);
	NEW_TMP_FILE(temp_filename)
	encoder.storeLibSVMProblem(temp_filename, problem);
	TEST_FILE_EQUAL(OPENMS_GET_TEST_DATA_PATH("LibSVMEncoder_test.txt"), temp_filename.c_str())

END_SECTION

START_SECTION((svm_problem* loadLibSVMProblem(const String& filename)))
	String allowed_characters = "ACNGT";
	svm_problem* problem;
	String temp_filename = OPENMS_GET_TEST_DATA_PATH("LibSVMEncoder_test.tmp");
	
	NEW_TMP_FILE(temp_filename)
	problem = encoder.loadLibSVMProblem(OPENMS_GET_TEST_DATA_PATH("LibSVMEncoder_test.txt"));
	encoder.storeLibSVMProblem(temp_filename, problem);
	TEST_FILE_EQUAL(OPENMS_GET_TEST_DATA_PATH("LibSVMEncoder_test.txt"), temp_filename.c_str())

END_SECTION

START_SECTION((void encodeOligoBorders(String sequence, UInt k_mer_length, const String& allowed_characters, UInt border_length, std::vector< std::pair<Int, DoubleReal> >& libsvm_vector, bool strict = false, bool unpaired=false, bool length_encoding = false)))
	String sequence = "ACNNGTATCA";
	String allowed_characters = "ACNGT";
	String output;
	UInt border_length = 3;
	vector< pair<Int, DoubleReal> > encoded_sequence;
	
	encoder.encodeOligoBorders(sequence, 1, allowed_characters, border_length, encoded_sequence);
	encoder.libSVMVectorToString(encoder.encodeLibSVMVector(encoded_sequence), output);
	TEST_EQUAL(output, "(2, 1) (2, 1) (3, 2) (3, 2) (4, 3) (6, 3) ")
	encoder.encodeOligoBorders(sequence, 2, allowed_characters, border_length, encoded_sequence);
	encoder.libSVMVectorToString(encoder.encodeLibSVMVector(encoded_sequence), output);
	TEST_EQUAL(output, "(3, 1) (3, 1) (9, 2) (11, 2) (14, 3) (22, 3) ")
	sequence = "ACNNGTZTCA";
	encoder.encodeOligoBorders(sequence, 1, allowed_characters, border_length, encoded_sequence);
	TEST_EQUAL(encoded_sequence.size(), 0)
	
END_SECTION

START_SECTION((svm_problem* encodeLibSVMProblemWithOligoBorderVectors(const std::vector< String > &sequences, std::vector< DoubleReal > &labels, UInt k_mer_length, const String &allowed_characters, UInt border_length, bool strict=false, bool unpaired=false, bool length_encoding=false)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	String output;
	UInt border_length = 3;
	vector< pair<Int, DoubleReal> > encoded_sequence;
  vector<DoubleReal> labels;
  struct svm_problem* data;

  labels.push_back(1);
  labels.push_back(2);
  sequences.push_back("ACNNGTATCA");
  sequences.push_back("AACNNGTACCA");
	data = encoder.encodeLibSVMProblemWithOligoBorderVectors(sequences, labels, 1, allowed_characters, border_length);
	encoder.libSVMVectorToString(data->x[0], output);
	TEST_EQUAL(output, "(2, 1) (2, 1) (3, 2) (3, 2) (4, 3) (6, 3) ")
	encoder.libSVMVectorToString(data->x[1], output);
	TEST_EQUAL(output, "(2, 1) (2, 1) (2, 2) (3, 2) (3, 3) (3, 3) ")
END_SECTION

START_SECTION((void encodeProblemWithOligoBorderVectors(const std::vector< AASequence > &sequences, UInt k_mer_length, const String &allowed_characters, UInt border_length, std::vector< std::vector< std::pair< Int, DoubleReal > > > &vectors)))
	vector<AASequence> sequences;
	String allowed_characters = "ACNGT";
	String output;
	UInt border_length = 3;
	vector< pair<Int, DoubleReal> > encoded_sequence;
  vector< vector< pair<Int, DoubleReal> > > encoded_sequences;
  
  sequences.push_back(AASequence("ACNNGTATCA"));
  sequences.push_back(AASequence("AACNNGTACCA"));
	encoder.encodeProblemWithOligoBorderVectors(sequences, 1, allowed_characters, border_length, encoded_sequences);
  TEST_EQUAL(encoded_sequences[0].size(), 6)
  TEST_EQUAL(encoded_sequences[0][0].first, 1)
  TEST_REAL_SIMILAR(encoded_sequences[0][0].second, 0.)
  TEST_EQUAL(encoded_sequences[0][1].first, 1)
  TEST_REAL_SIMILAR(encoded_sequences[0][1].second, 0.)
  TEST_EQUAL(encoded_sequences[0][2].first, 2)
  TEST_REAL_SIMILAR(encoded_sequences[0][2].second, 1.)
  TEST_EQUAL(encoded_sequences[0][3].first, 2)
  TEST_REAL_SIMILAR(encoded_sequences[0][3].second, 1.)
  TEST_EQUAL(encoded_sequences[0][4].first, 3)
  TEST_REAL_SIMILAR(encoded_sequences[0][4].second, 2.)
  TEST_EQUAL(encoded_sequences[0][5].first, 3)
  TEST_REAL_SIMILAR(encoded_sequences[0][5].second, 4.)

  TEST_EQUAL(encoded_sequences[1][0].first, 1)
  TEST_REAL_SIMILAR(encoded_sequences[1][0].second, 0.)
  TEST_EQUAL(encoded_sequences[1][1].first, 1)
  TEST_REAL_SIMILAR(encoded_sequences[1][1].second, 0.)
  TEST_EQUAL(encoded_sequences[1][2].first, 2)
  TEST_REAL_SIMILAR(encoded_sequences[1][2].second, 0.)
  TEST_EQUAL(encoded_sequences[1][3].first, 2)
  TEST_REAL_SIMILAR(encoded_sequences[1][3].second, 1.)
  TEST_EQUAL(encoded_sequences[1][4].first, 3)
  TEST_REAL_SIMILAR(encoded_sequences[1][4].second, 1.)
  TEST_EQUAL(encoded_sequences[1][5].first, 3)
  TEST_REAL_SIMILAR(encoded_sequences[1][5].second, 1.)
  
END_SECTION

START_SECTION((void encodeOligo(const AASequence &sequence, UInt k_mer_length, const String &allowed_characters, std::vector< std::pair< Int, DoubleReal > > &values, bool is_right_border=false)))
	AASequence sequence = AASequence("ACNNGTATCA");
	String allowed_characters = "ACNGT";
	String output;
	vector< pair<Int, DoubleReal> > encoded_sequence;
  ModificationsDB* modifications = ModificationsDB::getInstance();
  bool right_border = true;
	
	encoder.encodeOligo(sequence, 1, allowed_characters, encoded_sequence);
  TEST_EQUAL(encoded_sequence[0].first, 1)
  TEST_REAL_SIMILAR(encoded_sequence[0].second, 0.)
  TEST_EQUAL(encoded_sequence[1].first, 7)
  TEST_REAL_SIMILAR(encoded_sequence[1].second, 0.)
  TEST_EQUAL(encoded_sequence[2].first, 10)
  TEST_REAL_SIMILAR(encoded_sequence[2].second, 0.)
  TEST_EQUAL(encoded_sequence[3].first, 2)
  TEST_REAL_SIMILAR(encoded_sequence[3].second, 1.)
  TEST_EQUAL(encoded_sequence[4].first, 9)
  TEST_REAL_SIMILAR(encoded_sequence[4].second, 1.)
  TEST_EQUAL(encoded_sequence[5].first, 3)
  TEST_REAL_SIMILAR(encoded_sequence[5].second, 2.)
  TEST_EQUAL(encoded_sequence[6].first, 4)
  TEST_REAL_SIMILAR(encoded_sequence[6].second, 2.)
  TEST_EQUAL(encoded_sequence[7].first, 5)
  TEST_REAL_SIMILAR(encoded_sequence[7].second, 3.)
  TEST_EQUAL(encoded_sequence[8].first, 6)
  TEST_REAL_SIMILAR(encoded_sequence[8].second, 4.)
  TEST_EQUAL(encoded_sequence[9].first, 8)
  TEST_REAL_SIMILAR(encoded_sequence[9].second, 4.)
          
	encoder.encodeOligo(sequence, 1, allowed_characters, encoded_sequence, right_border);
  TEST_EQUAL(encoded_sequence[0].first, 1)
  TEST_REAL_SIMILAR(encoded_sequence[0].second, 0.)
  TEST_EQUAL(encoded_sequence[1].first, 4)
  TEST_REAL_SIMILAR(encoded_sequence[1].second, 0.)
  TEST_EQUAL(encoded_sequence[2].first, 10)
  TEST_REAL_SIMILAR(encoded_sequence[2].second, 0.)
  TEST_EQUAL(encoded_sequence[3].first, 2)
  TEST_REAL_SIMILAR(encoded_sequence[3].second, 1.)
  TEST_EQUAL(encoded_sequence[4].first, 9)
  TEST_REAL_SIMILAR(encoded_sequence[4].second, 1.)
  TEST_EQUAL(encoded_sequence[5].first, 7)
  TEST_REAL_SIMILAR(encoded_sequence[5].second, 2.)
  TEST_EQUAL(encoded_sequence[6].first, 8)
  TEST_REAL_SIMILAR(encoded_sequence[6].second, 2.)
  TEST_EQUAL(encoded_sequence[7].first, 6)
  TEST_REAL_SIMILAR(encoded_sequence[7].second, 3.)
  TEST_EQUAL(encoded_sequence[8].first, 3)
  TEST_REAL_SIMILAR(encoded_sequence[8].second, 4.)
  TEST_EQUAL(encoded_sequence[9].first, 5)
  TEST_REAL_SIMILAR(encoded_sequence[9].second, 4.)
  
  sequence = AASequence("ACNN");       
	encoder.encodeOligo(sequence, 2, allowed_characters, encoded_sequence);
  TEST_EQUAL(encoded_sequence[0].first, 1)
  TEST_REAL_SIMILAR(encoded_sequence[0].second, 1.)
  TEST_EQUAL(encoded_sequence[1].first, 2)
  TEST_REAL_SIMILAR(encoded_sequence[1].second, allowed_characters.size() * (modifications->getNumberOfModifications() + 1) + 2.)
  TEST_EQUAL(encoded_sequence[2].first, 3)
  TEST_REAL_SIMILAR(encoded_sequence[2].second, 2 * allowed_characters.size() * (modifications->getNumberOfModifications() + 1) + 2.)

  sequence = AASequence("ACNN");       
	encoder.encodeOligo(sequence, 2, allowed_characters, encoded_sequence, right_border);
  TEST_EQUAL(encoded_sequence[0].first, 3)
  TEST_REAL_SIMILAR(encoded_sequence[0].second, allowed_characters.size() * (modifications->getNumberOfModifications() + 1.))
  TEST_EQUAL(encoded_sequence[1].first, 2)
  TEST_REAL_SIMILAR(encoded_sequence[1].second, 2 * allowed_characters.size() * (modifications->getNumberOfModifications() + 1) + 1.)
  TEST_EQUAL(encoded_sequence[2].first, 1)
  TEST_REAL_SIMILAR(encoded_sequence[2].second, 2 * allowed_characters.size() * (modifications->getNumberOfModifications() + 1) + 2.)

END_SECTION

START_SECTION((void libSVMVectorToString(svm_node* vector, String& output)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	vector< pair<Int, DoubleReal> > encoded_sequence;
	svm_node* nodes;
	vector<svm_node*>::iterator it;
	String output;	
	String correct_output = "(1, 0.1) (2, 0.2) (4, 0.3) (5, 0.4) ";
	
	sequences.push_back(String("ACCGGGTTTT"));			
			
	encoder.encodeCompositionVector(sequences[0], encoded_sequence, allowed_characters);
	nodes = encoder.encodeLibSVMVector(encoded_sequence);
	
	encoder.libSVMVectorToString(nodes, output);
	
	TEST_EQUAL(output, correct_output)
	
END_SECTION

START_SECTION((void libSVMVectorsToString(svm_problem* vector, String& output)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	String output;	
	String correct_output = "(1, 0.1) (2, 0.2) (4, 0.3) (5, 0.4) \n(1, 0.5) (2, 0.5) \n";	
	vector<DoubleReal> labels;
	svm_problem* problem;
	
	labels.push_back(2.1);
	labels.push_back(1.3);
	sequences.push_back(String("ACCGGGTTTT"));			
	sequences.push_back(String("ACCA"));
		
	problem = encoder.encodeLibSVMProblemWithCompositionVectors(sequences, labels, allowed_characters);			
	encoder.libSVMVectorsToString(problem, output);
	TEST_EQUAL(output, correct_output)	
END_SECTION

START_SECTION(static void destroyProblem(svm_problem *problem))
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
