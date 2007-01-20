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
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <svm.h>

///////////////////////////

#include <string>
#include <vector>

///////////////////////////

START_TEST(SVMWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SVMWrapper* ptr;
SVMWrapper svm;

CHECK(SVMWrapper())
	ptr = new SVMWrapper();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(double getDoubleParameter(SVM_parameter_type type))
	svm.setParameter(C, 1.0043);
	svm.setParameter(NU, 0.0523);
	svm.setParameter(P, 1.2319);
	
	TEST_REAL_EQUAL(svm.getDoubleParameter(C), 1.0043)
	TEST_REAL_EQUAL(svm.getDoubleParameter(NU), 0.0523)
	TEST_REAL_EQUAL(svm.getDoubleParameter(P), 1.2319)
RESULT

CHECK(double getSVRProbability())
	LibSVMEncoder encoder;
	vector< vector< pair<SignedInt, DoubleReal> > > vectors;		
	vector< pair<SignedInt, DoubleReal> > temp_vector;
	vector<svm_node*>* encoded_vectors;
	UnsignedInt count = 100;
	vector<DoubleReal> labels;
	svm_problem* problem;
	
	for(UnsignedInt j = 0; j < count; j++)
	{	
		temp_vector.clear();
		for(UnsignedInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((DoubleReal) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoded_vectors = encoder.encodeLibSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(*encoded_vectors, &labels);
	svm.setParameter(PROBABILITY, 1);
	svm.train(problem);
	TEST_EQUAL(svm.getSVRProbability() == 0, false)
RESULT

CHECK(int getIntParameter(SVM_parameter_type type))
	svm.setParameter(SVM_TYPE, EPSILON_SVR);
	svm.setParameter(KERNEL_TYPE, LINEAR);
	svm.setParameter(DEGREE, 2);

	TEST_EQUAL(svm.getIntParameter(SVM_TYPE), EPSILON_SVR);
	TEST_EQUAL(svm.getIntParameter(KERNEL_TYPE), LINEAR);
	TEST_EQUAL(svm.getIntParameter(DEGREE), 2);
RESULT

CHECK(int train(struct svm_problem* problem))
	svm_problem* problem = new svm_problem();
	UnsignedInt count = 4;
	svm_node** nodes = new svm_node*[count];
	DoubleReal* labels = new DoubleReal[count];
	
	for(UnsignedInt i = 0; i < count; i++)
	{
		nodes[i] = new svm_node[count];
		nodes[i][count - 1].index = -1;
		labels[i] = i * 2 / 3 + 0.03;
	}
	problem->x = nodes;
	problem->l = count;
	problem->y = labels;
	TEST_EQUAL(svm.train(problem), 1)
RESULT

CHECK(static std::vector<DoubleReal>* getLabels(svm_problem* problem))
	svm_problem* problem = new svm_problem();
	UnsignedInt count = 4;
	svm_node** nodes = new svm_node*[count];
	DoubleReal* labels = new DoubleReal[count];
	std::vector<DoubleReal>* label_vector1 = new vector<DoubleReal>();
	std::vector<DoubleReal>* label_vector2 = new vector<DoubleReal>();
	
	for(UnsignedInt i = 0; i < count; i++)
	{
		nodes[i] = new svm_node[count];
		nodes[i][count - 1].index = -1;
		labels[i] = i * 2 / 3 + 0.03;
		label_vector1->push_back(labels[i]);
	}
	problem->x = nodes;
	problem->l = count;
	problem->y = labels;
	
	label_vector2 = SVMWrapper::getLabels(problem);
	TEST_EQUAL(label_vector1->size(), label_vector2->size())
	for(UnsignedInt i = 0; i < label_vector2->size(); i++)
	{
		TEST_REAL_EQUAL((*label_vector1)[i], (*label_vector2)[i])
	}	
	delete label_vector1;
	delete label_vector2;
	delete problem;
RESULT

CHECK((static std::vector<svm_problem*>* createRandomPartitions(svm_problem* problem, UnsignedInt number)))
	 svm_problem* problem = new svm_problem();
	UnsignedInt count = 4;
	 svm_node** nodes = new svm_node*[count];
	DoubleReal* labels = new DoubleReal[count];
	std::vector<svm_problem*>* partitions;
		
	for(UnsignedInt i = 0; i < count; i++)
	{
		nodes[i] = new svm_node[count];
		nodes[i][count - 1].index = -1;
		labels[i] = i * 2 / 3 + 0.03;
	}
	problem->x = nodes;
	problem->l = count;
	problem->y = labels;
	
	partitions = SVMWrapper::createRandomPartitions(problem, 2);
	TEST_EQUAL(partitions->size(), 2)
	TEST_EQUAL((*partitions)[0]->l, 2)
	TEST_EQUAL((*partitions)[1]->l, 2)
RESULT

CHECK((static svm_problem* mergePartitions(const std::vector<svm_problem*>* const problems, UnsignedInt except)))
	 svm_problem* problem = new svm_problem();
	 svm_problem* problem2;
	 UnsignedInt count = 10;
	 UnsignedInt number_of_partitions = 5;
	 svm_node** nodes = new svm_node*[count];
	 DoubleReal* labels = new DoubleReal[count];
	 std::vector<svm_problem*>* partitions;
	
		
	for(UnsignedInt i = 0; i < count; i++)
	{
		nodes[i] = new svm_node[count];
		nodes[i][count - 1].index = -1;
		labels[i] = ((DoubleReal) i * 2) / 3 + 0.03;
		for(UnsignedInt j = 0; j < count; j++)
		{
			nodes[i][j].value = ((DoubleReal) i * 2) / 3;
		}
	}
	problem->x = nodes;
	problem->l = count;
	problem->y = labels;

	partitions = SVMWrapper::createRandomPartitions(problem, number_of_partitions);
	problem2 = SVMWrapper::mergePartitions(partitions, 4);
	UnsignedInt problem2_size = (count / number_of_partitions) * (number_of_partitions - 1);
	UnsignedInt partition_size = count / number_of_partitions;
	TEST_EQUAL((UnsignedInt) problem2->l, problem2_size)
	for(UnsignedInt i = 0; i < problem2_size; i++)
	{
		UnsignedInt j = 0;
		while(problem->x[i][j].index != -1 && problem2->x[i][j].index != -1)
		{
			TEST_REAL_EQUAL((*partitions)[i / partition_size]->x[i % partition_size][j].value, problem2->x[i][j].value)
			++j;
		}
		TEST_REAL_EQUAL((*partitions)[i / partition_size]->y[i % partition_size], problem2->y[i])
	}
RESULT

CHECK((static void calculateGaussTable(UnsignedInt border_length, DoubleReal sigma, std::vector<DoubleReal>& gauss_table)))
  UnsignedInt border_length = 5;
  DoubleReal sigma = 2;
  DoubleReal sigma_square = sigma * sigma;
  vector<DoubleReal> gauss_table;
  svm.calculateGaussTable(border_length, sigma, gauss_table);
  
  TEST_EQUAL(gauss_table.size(), 5)
  TEST_EQUAL(gauss_table[0], 1)
  TEST_EQUAL(gauss_table[1], exp((-1 / (4.0 * sigma_square)) * 1))
  TEST_EQUAL(gauss_table[2], exp((-1 / (4.0 * sigma_square)) * 4))
  TEST_EQUAL(gauss_table[3], exp((-1 / (4.0 * sigma_square)) * 9))
  TEST_EQUAL(gauss_table[4], exp((-1 / (4.0 * sigma_square)) * 16))	
RESULT

CHECK((std::map<SVM_parameter_type, DoubleReal>* performCrossValidation(svm_problem* problem, std::map<SVM_parameter_type, DoubleReal>& start_values, std::map<SVM_parameter_type, DoubleReal>& step_sizes, std::map<SVM_parameter_type, DoubleReal>& end_values, DoubleReal* cv_quality, UnsignedInt number_of_partitions, UnsignedInt number_of_runs, bool additive_step_size = true, bool output = false, String performances_file_name = "performances.txt")))
	map<SVM_parameter_type, DoubleReal> start_values;
	map<SVM_parameter_type, DoubleReal> step_sizes;
	map<SVM_parameter_type, DoubleReal> end_values;
	LibSVMEncoder encoder;
	vector< vector< pair<SignedInt, DoubleReal> > > vectors;		
	vector< pair<SignedInt, DoubleReal> > temp_vector;
	vector<svm_node*>* encoded_vectors;
	UnsignedInt count = 8;
	vector<DoubleReal> labels;
	svm_problem* problem;
	map<SVM_parameter_type, DoubleReal>* parameters;
	DoubleReal cv_quality;
	
	for(UnsignedInt j = 0; j < count; j++)
	{	
		temp_vector.clear();
		for(UnsignedInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((DoubleReal) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoded_vectors = encoder.encodeLibSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(*encoded_vectors, &labels);

	start_values.insert(make_pair(C, 1));
	step_sizes.insert(make_pair(C, 100));
	end_values.insert(make_pair(C, 1000));

	start_values.insert(make_pair(NU, 0.4));
	step_sizes.insert(make_pair(NU, 0.1));
	end_values.insert(make_pair(NU, 0.6));

	start_values.insert(make_pair(DEGREE, 1));
	step_sizes.insert(make_pair(DEGREE, 1));
	end_values.insert(make_pair(DEGREE, 3));

	parameters = svm.performCrossValidation(problem, start_values, step_sizes, end_values, &cv_quality, 2, 1, true, false);
	TEST_NOT_EQUAL(parameters->size(), 0)
RESULT

CHECK(std::vector<DoubleReal>* predict(struct svm_problem* predictProblem))
 	LibSVMEncoder encoder;
	vector< vector< pair<SignedInt, DoubleReal> > > vectors;		
	vector< pair<SignedInt, DoubleReal> > temp_vector;
	vector<svm_node*>* encoded_vectors;
	UnsignedInt count = 8;
	vector<DoubleReal> labels;
	vector<DoubleReal>* predicted_labels;
	svm_problem* problem;
	
	for(UnsignedInt j = 0; j < count; j++)
	{	
		temp_vector.clear();
		for(UnsignedInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((DoubleReal) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoded_vectors = encoder.encodeLibSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(*encoded_vectors, &labels);
	svm.train(problem);
	predicted_labels = svm.predict(problem);
	TEST_NOT_EQUAL(predicted_labels->size(), 0)
RESULT

CHECK((svm_problem* computeKernelMatrix(svm_problem* problem1, svm_problem* problem2)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	String output;
	SignedInt border_length = 5;
	DoubleReal sigma = 2;
	vector< pair<SignedInt, DoubleReal> > encoded_sequence;
  vector<DoubleReal> labels;
  struct svm_problem* data;
  struct svm_problem* kernel_matrix;
  LibSVMEncoder encoder;
	
	svm.setParameter(BORDER_LENGTH, border_length);
	svm.setParameter(SIGMA, sigma);
	svm.setParameter(KERNEL_TYPE, OLIGO);
  labels.push_back(1);
  labels.push_back(2);
  sequences.push_back("ACNNGTATCA");
  sequences.push_back("AACNNGTACCA");
	data = encoder.encodeLibSVMProblemWithOligoBorderVectors(sequences, &labels, 1, allowed_characters, border_length);
	kernel_matrix = svm.computeKernelMatrix(data, data);
	svm.train(data);

	PRECISION(0.0001)
	TEST_REAL_EQUAL(kernel_matrix->x[0][0].value, 1)
	TEST_REAL_EQUAL(kernel_matrix->x[0][1].value, 19.7156)
	TEST_REAL_EQUAL(kernel_matrix->x[0][2].value, 21.1308)
	TEST_REAL_EQUAL(kernel_matrix->x[1][0].value, 2)
	TEST_REAL_EQUAL(kernel_matrix->x[1][1].value, 21.1308)
	TEST_REAL_EQUAL(kernel_matrix->x[1][2].value, 27.2309)
	TEST_EQUAL(kernel_matrix->x[0][0].index, 0)
	TEST_EQUAL(kernel_matrix->x[0][1].index, 1)
	TEST_EQUAL(kernel_matrix->x[0][2].index, 2)
	TEST_EQUAL(kernel_matrix->x[1][0].index, 0)
	TEST_EQUAL(kernel_matrix->x[1][1].index, 1)
	TEST_EQUAL(kernel_matrix->x[1][2].index, 2)
	TEST_EQUAL(kernel_matrix->y[0], 1)
	TEST_EQUAL(kernel_matrix->y[1], 2)

RESULT

CHECK((static DoubleReal kernelOligo(const svm_node* x, const svm_node* y, const std::vector<DoubleReal>& gauss_table, DoubleReal sigma_square = 0, UnsignedInt max_distance = 50)))
  vector<DoubleReal> labels;
	String sequence = "ACNNGTATCA";
	String allowed_characters = "ACNGT";
	String output;
	SignedInt border_length = 5;
	svm_problem* data;
	DoubleReal result = 0;
  DoubleReal sigma = 2;
  vector<DoubleReal> gauss_table;
	vector<String> sequences;
  svm.calculateGaussTable(border_length, sigma, gauss_table);
	LibSVMEncoder encoder;
	svm.setParameter(BORDER_LENGTH, border_length);
	svm.setParameter(SIGMA, sigma);
	svm.setParameter(KERNEL_TYPE, OLIGO);
	
  labels.push_back(1);
  labels.push_back(2);
  sequences.push_back("ACNNGTATCA");
  sequences.push_back("AACNNGTACCA");
	data = encoder.encodeLibSVMProblemWithOligoBorderVectors(sequences, &labels, 1, allowed_characters, border_length);
	result = SVMWrapper::kernelOligo(data->x[0], data->x[1], gauss_table);
	PRECISION(0.0001)
	TEST_REAL_EQUAL(result, 21.1308)
	delete data;
RESULT

CHECK((void getDecisionValues(svm_problem* data, std::vector<DoubleReal>& decision_values)))
 	LibSVMEncoder encoder;
	vector< vector< pair<SignedInt, DoubleReal> > > vectors;		
	vector< pair<SignedInt, DoubleReal> > temp_vector;
	vector<svm_node*>* encoded_vectors;
	UnsignedInt count = 8;
	vector<DoubleReal> labels;
	vector<DoubleReal>* predicted_labels;
	svm_problem* problem;
	vector<DoubleReal> decision_values;
	
	svm.setParameter(SVM_TYPE, NU_SVR);
	svm.setParameter(KERNEL_TYPE, POLY);
	svm.setParameter(DEGREE, 2);
	for(UnsignedInt j = 0; j < count; j++)
	{	
		temp_vector.clear();
		for(UnsignedInt i = 1; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((DoubleReal) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoded_vectors = encoder.encodeLibSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(*encoded_vectors, &labels);
	svm.train(problem);
	predicted_labels = svm.predict(problem);
	TEST_NOT_EQUAL(predicted_labels->size(), 0)
	svm.getDecisionValues(problem, decision_values);
	TEST_EQUAL(*predicted_labels == decision_values, true)
	
RESULT

CHECK((void scaleData(svm_problem* data, SignedInt max_scale_value = -1)))
 	LibSVMEncoder encoder;
	vector< vector< pair<SignedInt, DoubleReal> > > vectors;		
	vector< pair<SignedInt, DoubleReal> > temp_vector;
	vector<svm_node*>* encoded_vectors;
	UnsignedInt count = 8;
	vector<DoubleReal> labels;
	svm_problem* problem;
	vector<DoubleReal> decision_values;
	
	svm.setParameter(SVM_TYPE, NU_SVR);
	svm.setParameter(KERNEL_TYPE, POLY);
	svm.setParameter(DEGREE, 2);
	for(UnsignedInt j = 0; j < count; j++)
	{	
		temp_vector.clear();
		for(UnsignedInt i = 1; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i, ((DoubleReal) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoded_vectors = encoder.encodeLibSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(*encoded_vectors, &labels);
	svm.scaleData(problem, 2);

	TEST_REAL_EQUAL(problem->x[0][0].value, 0)
	TEST_REAL_EQUAL(problem->x[0][1].value, 0)
	TEST_REAL_EQUAL(problem->x[0][2].value, 0)
	TEST_REAL_EQUAL(problem->x[0][3].value, 0)
	TEST_REAL_EQUAL(problem->x[0][4].value, 0)
	TEST_REAL_EQUAL(problem->x[1][0].value, 0.2857)
	TEST_REAL_EQUAL(problem->x[1][1].value, 0.2857)
	TEST_REAL_EQUAL(problem->x[1][2].value, 0.2857)
	TEST_REAL_EQUAL(problem->x[1][3].value, 0.2857)
	TEST_REAL_EQUAL(problem->x[1][4].value, 0.2857)
	TEST_REAL_EQUAL(problem->x[2][0].value, 0.5714)
	TEST_REAL_EQUAL(problem->x[2][1].value, 0.5714)
	TEST_REAL_EQUAL(problem->x[2][2].value, 0.5714)
	TEST_REAL_EQUAL(problem->x[2][3].value, 0.5714)
	TEST_REAL_EQUAL(problem->x[2][4].value, 0.5714)
	TEST_REAL_EQUAL(problem->x[3][0].value, 0.8571)
	TEST_REAL_EQUAL(problem->x[3][1].value, 0.8571)
	TEST_REAL_EQUAL(problem->x[3][2].value, 0.8571)
	TEST_REAL_EQUAL(problem->x[3][3].value, 0.8571)
	TEST_REAL_EQUAL(problem->x[3][4].value, 0.8571)
	TEST_REAL_EQUAL(problem->x[4][0].value, 1.1429)
	TEST_REAL_EQUAL(problem->x[4][1].value, 1.1429)
	TEST_REAL_EQUAL(problem->x[4][2].value, 1.1429)
	TEST_REAL_EQUAL(problem->x[4][3].value, 1.1429)
	TEST_REAL_EQUAL(problem->x[4][4].value, 1.1429)
	TEST_REAL_EQUAL(problem->x[5][0].value, 1.4286)
	TEST_REAL_EQUAL(problem->x[5][1].value, 1.4286)
	TEST_REAL_EQUAL(problem->x[5][2].value, 1.4286)
	TEST_REAL_EQUAL(problem->x[5][3].value, 1.4286)
	TEST_REAL_EQUAL(problem->x[5][4].value, 1.4286)
	TEST_REAL_EQUAL(problem->x[6][0].value, 1.7143)
	TEST_REAL_EQUAL(problem->x[6][1].value, 1.7143)
	TEST_REAL_EQUAL(problem->x[6][2].value, 1.7143)
	TEST_REAL_EQUAL(problem->x[6][3].value, 1.7143)
	TEST_REAL_EQUAL(problem->x[6][4].value, 1.7143)
	TEST_REAL_EQUAL(problem->x[7][0].value, 2)
	TEST_REAL_EQUAL(problem->x[7][1].value, 2)
	TEST_REAL_EQUAL(problem->x[7][2].value, 2)
	TEST_REAL_EQUAL(problem->x[7][3].value, 2)
	TEST_REAL_EQUAL(problem->x[7][4].value, 2)

	svm.scaleData(problem);

	TEST_REAL_EQUAL(problem->x[0][0].value, -1)
	TEST_REAL_EQUAL(problem->x[0][1].value, -1)
	TEST_REAL_EQUAL(problem->x[0][2].value, -1)
	TEST_REAL_EQUAL(problem->x[0][3].value, -1)
	TEST_REAL_EQUAL(problem->x[0][4].value, -1)
	TEST_REAL_EQUAL(problem->x[1][0].value, -0.7143)
	TEST_REAL_EQUAL(problem->x[1][1].value, -0.7143)
	TEST_REAL_EQUAL(problem->x[1][2].value, -0.7143)
	TEST_REAL_EQUAL(problem->x[1][3].value, -0.7143)
	TEST_REAL_EQUAL(problem->x[1][4].value, -0.7143)
	TEST_REAL_EQUAL(problem->x[2][0].value, -0.4286)
	TEST_REAL_EQUAL(problem->x[2][1].value, -0.4286)
	TEST_REAL_EQUAL(problem->x[2][2].value, -0.4286)
	TEST_REAL_EQUAL(problem->x[2][3].value, -0.4286)
	TEST_REAL_EQUAL(problem->x[2][4].value, -0.4286)
	TEST_REAL_EQUAL(problem->x[3][0].value, -0.1429)
	TEST_REAL_EQUAL(problem->x[3][1].value, -0.1429)
	TEST_REAL_EQUAL(problem->x[3][2].value, -0.1429)
	TEST_REAL_EQUAL(problem->x[3][3].value, -0.1429)
	TEST_REAL_EQUAL(problem->x[3][4].value, -0.1429)
	TEST_REAL_EQUAL(problem->x[4][0].value, 0.1429)
	TEST_REAL_EQUAL(problem->x[4][1].value, 0.1429)
	TEST_REAL_EQUAL(problem->x[4][2].value, 0.1429)
	TEST_REAL_EQUAL(problem->x[4][3].value, 0.1429)
	TEST_REAL_EQUAL(problem->x[4][4].value, 0.1429)
	TEST_REAL_EQUAL(problem->x[5][0].value, 0.4286)
	TEST_REAL_EQUAL(problem->x[5][1].value, 0.4286)
	TEST_REAL_EQUAL(problem->x[5][2].value, 0.4286)
	TEST_REAL_EQUAL(problem->x[5][3].value, 0.4286)
	TEST_REAL_EQUAL(problem->x[5][4].value, 0.4286)
	TEST_REAL_EQUAL(problem->x[6][0].value, 0.7143)
	TEST_REAL_EQUAL(problem->x[6][1].value, 0.7143)
	TEST_REAL_EQUAL(problem->x[6][2].value, 0.7143)
	TEST_REAL_EQUAL(problem->x[6][3].value, 0.7143)
	TEST_REAL_EQUAL(problem->x[6][4].value, 0.7143)
	TEST_REAL_EQUAL(problem->x[7][0].value, 1)
	TEST_REAL_EQUAL(problem->x[7][1].value, 1)
	TEST_REAL_EQUAL(problem->x[7][2].value, 1)
	TEST_REAL_EQUAL(problem->x[7][3].value, 1)
	TEST_REAL_EQUAL(problem->x[7][4].value, 1)

RESULT

CHECK((void getSignificanceBorders(svm_problem* data, std::pair<DoubleReal, DoubleReal>& borders, DoubleReal confidence = 0.95, UnsignedInt number_of_runs = 10, UnsignedInt number_of_partitions = 5, DoubleReal step_size = 0.01, UnsignedInt max_iterations = 1000000)))
  // ???
RESULT

CHECK((DoubleReal getPValue(DoubleReal sigma1, DoubleReal sigma2, std::pair<DoubleReal, DoubleReal> point)))
  // ???
RESULT

CHECK(void setTrainingSample(svm_problem* training_sample))
  // ???
RESULT

CHECK((void setParameter(SVM_parameter_type type, double value)))
 	svm.setParameter(C, 1.0043);
	svm.setParameter(NU, 0.0523);
	svm.setParameter(P, 1.2319);
	
	TEST_REAL_EQUAL(svm.getDoubleParameter(C), 1.0043)
	TEST_REAL_EQUAL(svm.getDoubleParameter(NU), 0.0523)
	TEST_REAL_EQUAL(svm.getDoubleParameter(P), 1.2319)
RESULT

CHECK((void setParameter(SVM_parameter_type type, int value)))
	svm.setParameter(SVM_TYPE, EPSILON_SVR);
	svm.setParameter(KERNEL_TYPE, LINEAR);
	svm.setParameter(DEGREE, 2);
	svm.setParameter(C, 23);
	svm.setParameter(PROBABILITY, 1);

	TEST_EQUAL(svm.getIntParameter(SVM_TYPE), EPSILON_SVR);
	TEST_EQUAL(svm.getIntParameter(KERNEL_TYPE), LINEAR);
	TEST_EQUAL(svm.getIntParameter(DEGREE), 2);
	TEST_EQUAL((int) svm.getDoubleParameter(C), 23);
	TEST_EQUAL(svm.getIntParameter(PROBABILITY), 1)
RESULT

CHECK(~SVMWrapper())
	delete ptr;
RESULT

CHECK((void loadModel(std::string modelFilename)))
	LibSVMEncoder encoder;
	svm.setParameter(KERNEL_TYPE, POLY);
	vector< vector< pair<SignedInt, DoubleReal> > > vectors;		
	vector< pair<SignedInt, DoubleReal> > temp_vector;
	vector<svm_node*>* encoded_vectors;
	UnsignedInt count = 8;
	vector<DoubleReal> labels;
	vector<DoubleReal>* predicted_labels1;
	vector<DoubleReal>* predicted_labels2;
	svm_problem* problem;
	SVMWrapper svm2;
	
	for(UnsignedInt j = 0; j < count; j++)
	{	
		temp_vector.clear();
		for(UnsignedInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((DoubleReal) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoded_vectors = encoder.encodeLibSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(*encoded_vectors, &labels);
	svm.train(problem);
	predicted_labels1 = svm.predict(problem);
	
	String filename = "svm.model";
	NEW_TMP_FILE(filename)
	svm.saveModel(filename);
	svm2.loadModel(filename);
	predicted_labels2 = svm2.predict(problem);
	TEST_NOT_EQUAL(predicted_labels1->size(), 0)
	TEST_EQUAL(predicted_labels1->size(), predicted_labels2->size())
	for(UnsignedInt i = 0; i < predicted_labels1->size(); i++)
	{
		TEST_REAL_EQUAL((*predicted_labels1)[i], (*predicted_labels2)[i])
	}
RESULT

CHECK((void saveModel(std::string modelFilename) const throw(Exception::UnableToCreateFile)))
	LibSVMEncoder encoder;
	svm.setParameter(KERNEL_TYPE, POLY);
	vector< vector< pair<SignedInt, DoubleReal> > > vectors;		
	vector< pair<SignedInt, DoubleReal> > temp_vector;
	vector<svm_node*>* encoded_vectors;
	UnsignedInt count = 8;
	vector<DoubleReal> labels;
	vector<DoubleReal>* predicted_labels1;
	vector<DoubleReal>* predicted_labels2;
	svm_problem* problem;
	SVMWrapper svm2;	
	
	for(UnsignedInt j = 0; j < count; j++)
	{	
		temp_vector.clear();
		for(UnsignedInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((DoubleReal) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoded_vectors = encoder.encodeLibSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(*encoded_vectors, &labels);
	svm.train(problem);
	
	String filename = "svm.model";
	NEW_TMP_FILE(filename)
	svm.saveModel(filename);
	svm2.loadModel(filename);
	predicted_labels1 = svm.predict(problem);
	predicted_labels2 = svm2.predict(problem);
	TEST_NOT_EQUAL(predicted_labels1->size(), 0)
	TEST_NOT_EQUAL(predicted_labels2->size(), 0)
	TEST_EQUAL(predicted_labels1->size(), predicted_labels2->size())

	for(UnsignedInt i = 0; i < predicted_labels1->size(); i++)
	{
		TEST_REAL_EQUAL((*predicted_labels1)[i], (*predicted_labels2)[i])
	}
RESULT

CHECK((std::vector<DoubleReal>* predict(const std::vector<svm_node*>& vectors)))
	LibSVMEncoder encoder;
	vector< vector< pair<SignedInt, DoubleReal> > > vectors;		
	vector< pair<SignedInt, DoubleReal> > temp_vector;
	vector<svm_node*>* encoded_vectors;
	UnsignedInt count = 8;
	vector<DoubleReal> labels;
	vector<DoubleReal>* predicted_labels;
	svm_problem* problem;
	
	for(UnsignedInt j = 0; j < count; j++)
	{	
		for(UnsignedInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((DoubleReal) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoded_vectors = encoder.encodeLibSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(*encoded_vectors, &labels);
	svm.train(problem);
	predicted_labels = svm.predict(*encoded_vectors);
	TEST_NOT_EQUAL(predicted_labels->size(), 0)
	
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
