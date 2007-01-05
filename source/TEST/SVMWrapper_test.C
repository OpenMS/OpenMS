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

CHECK((SVMWrapper()))
	ptr = new SVMWrapper();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~SVMWrapper()))
	delete ptr;
RESULT

CHECK((double getDoubleParameter(SVM_parameter_type type)))
	svm.setParameter(C, 1.0043);
	svm.setParameter(NU, 0.0523);
	svm.setParameter(P, 1.2319);
	
	TEST_REAL_EQUAL(svm.getDoubleParameter(C), 1.0043)
	TEST_REAL_EQUAL(svm.getDoubleParameter(NU), 0.0523)
	TEST_REAL_EQUAL(svm.getDoubleParameter(P), 1.2319)
	
RESULT

CHECK((void setParameter(SVM_parameter_type type, double value)))
	svm.setParameter(C, 1.0043);
	svm.setParameter(NU, 0.0523);
	svm.setParameter(P, 1.2319);
	
	TEST_REAL_EQUAL(svm.getDoubleParameter(C), 1.0043)
	TEST_REAL_EQUAL(svm.getDoubleParameter(NU), 0.0523)
	TEST_REAL_EQUAL(svm.getDoubleParameter(P), 1.2319)
RESULT

CHECK((int getIntParameter(SVM_parameter_type type)))
	svm.setParameter(SVM_TYPE, EPSILON_SVR);
	svm.setParameter(KERNEL_TYPE, LINEAR);
	svm.setParameter(DEGREE, 2);

	TEST_EQUAL(svm.getIntParameter(SVM_TYPE), EPSILON_SVR);
	TEST_EQUAL(svm.getIntParameter(KERNEL_TYPE), LINEAR);
	TEST_EQUAL(svm.getIntParameter(DEGREE), 2);
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

CHECK((double getSVRProbability()))
	LibSVMEncoder encoder;
	vector< vector< pair<UnsignedInt, DoubleReal> > > vectors;		
	vector< pair<UnsignedInt, DoubleReal> > temp_vector;
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
	encoded_vectors = encoder.encodeLIBSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLIBSVMProblem(*encoded_vectors, &labels);
	svm.setParameter(PROBABILITY, 1);
	svm.train(problem);
	TEST_EQUAL(svm.getSVRProbability() == 0, false)
RESULT

CHECK((int train(struct svm_problem* problem)))
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

CHECK((static std::vector<DoubleReal>* getLabels(svm_problem* problem)))
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

CHECK((std::map<SVM_parameter_type, DoubleReal>* performCrossValidation(svm_problem* problem, std::map<SVM_parameter_type, DoubleReal>& start_values, std::map<SVM_parameter_type, DoubleReal>& step_sizes, std::map<SVM_parameter_type, DoubleReal>& end_values, DoubleReal* cv_quality, UnsignedInt number_of_partitions, UnsignedInt number_of_runs, bool additive_step_size = true, bool output = false)))
	map<SVM_parameter_type, DoubleReal> start_values;
	map<SVM_parameter_type, DoubleReal> step_sizes;
	map<SVM_parameter_type, DoubleReal> end_values;
	LibSVMEncoder encoder;
	vector< vector< pair<UnsignedInt, DoubleReal> > > vectors;		
	vector< pair<UnsignedInt, DoubleReal> > temp_vector;
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
	encoded_vectors = encoder.encodeLIBSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLIBSVMProblem(*encoded_vectors, &labels);

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

CHECK((std::vector<DoubleReal>* predict(const std::vector<svm_node*>& vectors)))
	LibSVMEncoder encoder;
	vector< vector< pair<UnsignedInt, DoubleReal> > > vectors;		
	vector< pair<UnsignedInt, DoubleReal> > temp_vector;
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
	encoded_vectors = encoder.encodeLIBSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLIBSVMProblem(*encoded_vectors, &labels);
	svm.train(problem);
	predicted_labels = svm.predict(*encoded_vectors);
	TEST_NOT_EQUAL(predicted_labels->size(), 0)
	
RESULT

CHECK((std::vector<DoubleReal>* predict(struct svm_problem* predictProblem)))
	LibSVMEncoder encoder;
	vector< vector< pair<UnsignedInt, DoubleReal> > > vectors;		
	vector< pair<UnsignedInt, DoubleReal> > temp_vector;
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
	encoded_vectors = encoder.encodeLIBSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLIBSVMProblem(*encoded_vectors, &labels);
	svm.train(problem);
	predicted_labels = svm.predict(problem);
	TEST_NOT_EQUAL(predicted_labels->size(), 0)

RESULT

CHECK((void loadModel(std::string modelFilename)))
	LibSVMEncoder encoder;
	vector< vector< pair<UnsignedInt, DoubleReal> > > vectors;		
	vector< pair<UnsignedInt, DoubleReal> > temp_vector;
	vector<svm_node*>* encoded_vectors;
	UnsignedInt count = 8;
	vector<DoubleReal> labels;
	vector<DoubleReal>* predicted_labels1;
	vector<DoubleReal>* predicted_labels2;
	svm_problem* problem;
	SVMWrapper svm2;

  svm.setParameter(KERNEL_TYPE, POLY);
  svm2.setParameter(KERNEL_TYPE, POLY);	
	
	for(UnsignedInt j = 0; j < count; j++)
	{	
		temp_vector.clear();
		for(UnsignedInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((DoubleReal) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoded_vectors = encoder.encodeLIBSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLIBSVMProblem(*encoded_vectors, &labels);
	svm.train(problem);
	predicted_labels1 = svm.predict(problem);
	
	String filename = "svm.model";
	NEW_TMP_FILE(filename)
	svm.saveModel(filename);
	svm2.loadModel(filename);
	TEST_EQUAL(svm.getIntParameter(KERNEL_TYPE), POLY);
	TEST_EQUAL(svm2.getIntParameter(KERNEL_TYPE), POLY);
	predicted_labels2 = svm2.predict(problem);
	TEST_NOT_EQUAL(predicted_labels1->size(), 0)
	TEST_EQUAL(predicted_labels1->size(), predicted_labels2->size())
	for(UnsignedInt i = 0; i < predicted_labels1->size(); i++)
	{
		TEST_REAL_EQUAL((*predicted_labels1)[i], (*predicted_labels2)[i])
	}
RESULT

CHECK((void saveModel(std::string modelFilename)))
	LibSVMEncoder encoder;
	vector< vector< pair<UnsignedInt, DoubleReal> > > vectors;		
	vector< pair<UnsignedInt, DoubleReal> > temp_vector;
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
	encoded_vectors = encoder.encodeLIBSVMVectors(vectors);
	for(UnsignedInt i = 0; i < count; i++)
	{
		labels.push_back(((DoubleReal) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLIBSVMProblem(*encoded_vectors, &labels);
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
