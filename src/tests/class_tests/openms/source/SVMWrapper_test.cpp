// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
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
SVMWrapper* nullPointer = nullptr;
SVMWrapper svm;

START_SECTION((SVMWrapper()))
	ptr = new SVMWrapper();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((double getDoubleParameter(SVM_parameter_type type)))
	svm.setParameter(SVMWrapper::C, 1.0043);
	svm.setParameter(SVMWrapper::NU, 0.0523);
	svm.setParameter(SVMWrapper::P, 1.2319);

	TEST_REAL_SIMILAR(svm.getDoubleParameter(SVMWrapper::C), 1.0043)
	TEST_REAL_SIMILAR(svm.getDoubleParameter(SVMWrapper::NU), 0.0523)
	TEST_REAL_SIMILAR(svm.getDoubleParameter(SVMWrapper::P), 1.2319)
END_SECTION

START_SECTION((double getSVRProbability()))
	LibSVMEncoder encoder;
	vector< vector< pair<Int, double> > > vectors;
	vector< pair<Int, double> > temp_vector;
	vector<svm_node*> encoded_vectors;
	Int count = 100;
	vector<double> labels;
	svm_problem* problem;

	for (Int j = 0; j < count; j++)
	{
		temp_vector.clear();
		for (Int i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((double) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoder.encodeLibSVMVectors(vectors, encoded_vectors);
	for (Int i = 0; i < count; i++)
	{
		labels.push_back(((double) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.setParameter(SVMWrapper::PROBABILITY, 1);
	svm.train(problem);
	TEST_EQUAL(svm.getSVRProbability() == 0, false)
END_SECTION

START_SECTION((Int getIntParameter(SVM_parameter_type type)))
	svm.setParameter(SVMWrapper::SVM_TYPE, EPSILON_SVR);
	svm.setParameter(SVMWrapper::KERNEL_TYPE, LINEAR);
	svm.setParameter(SVMWrapper::DEGREE, 2);

	TEST_EQUAL(svm.getIntParameter(SVMWrapper::SVM_TYPE)==EPSILON_SVR,true);
	TEST_EQUAL(svm.getIntParameter(SVMWrapper::KERNEL_TYPE)==LINEAR,true);
	TEST_EQUAL(svm.getIntParameter(SVMWrapper::DEGREE)==2,true);
END_SECTION

START_SECTION((Int train(struct svm_problem *problem)))
	svm_problem* problem = new svm_problem();
	UInt count = 4;
	svm_node** nodes = new svm_node*[count];
	double* labels = new double[count];

	for (Size i = 0; i < count; i++)
	{
		nodes[i] = new svm_node[count];
		nodes[i][count - 1].index = -1;
		labels[i] = i * 2 / 3 + 0.03;
	}
	problem->x = nodes;
	problem->l = count;
	problem->y = labels;
	TEST_EQUAL(svm.train(problem), 1)
END_SECTION

START_SECTION((Int train(SVMData &problem)))

	SVMWrapper svm2;
	SVMData problem;
	UInt count = 4;
	vector<double> labels;
	vector< vector<pair<Int, double> > > sequences;
	vector<pair<Int, double> > sequence;
	
	svm2.setParameter(SVMWrapper::KERNEL_TYPE, SVMWrapper::OLIGO);
	svm2.setParameter(SVMWrapper::BORDER_LENGTH, 2);
	svm2.setParameter(SVMWrapper::C, 1);
	svm2.setParameter(SVMWrapper::SIGMA, 1);
	svm2.setParameter(SVMWrapper::SVM_TYPE, NU_SVR);

	for (Size i = 0; i < count; i++)
	{
		sequence.clear();
		sequence.push_back(make_pair(1, rand()));
		sequence.push_back(make_pair(2, rand()));
		sequence.push_back(make_pair(3, rand()));
		sequence.push_back(make_pair(4, rand()));
		sequences.push_back(sequence);
		labels.push_back(i * 2 / 3 + 0.03);
	}
	problem.sequences = sequences;
	problem.labels = labels;

	TEST_EQUAL(svm2.train(problem), 1)
END_SECTION

START_SECTION((static void getLabels(svm_problem *problem, std::vector< double > &labels)))
	svm_problem* problem = new svm_problem();
	UInt count = 4;
	svm_node** nodes = new svm_node*[count];
	double* labels = new double[count];
	std::vector<double> label_vector1;
	std::vector<double> label_vector2;

	for (Size i = 0; i < count; i++)
	{
		nodes[i] = new svm_node[count];
		nodes[i][count - 1].index = -1;
		labels[i] = i * 2 / 3 + 0.03;
		label_vector1.push_back(labels[i]);
	}
	problem->x = nodes;
	problem->l = count;
	problem->y = labels;

	SVMWrapper::getLabels(problem, label_vector2);
	TEST_EQUAL(label_vector1.size(), label_vector2.size())
	for (Size i = 0; i < label_vector2.size(); i++)
	{
		TEST_REAL_SIMILAR(label_vector1[i], label_vector2[i])
	}
	delete problem;
END_SECTION

START_SECTION((static void createRandomPartitions(svm_problem *problem, Size number, std::vector< svm_problem * > &partitions)))
	 svm_problem* problem = new svm_problem();
	UInt count = 4;
	 svm_node** nodes = new svm_node*[count];
	double* labels = new double[count];
	std::vector<svm_problem*> partitions;

	for (Size i = 0; i < count; i++)
	{
		nodes[i] = new svm_node[count];
		nodes[i][count - 1].index = -1;
		labels[i] = i * 2 / 3 + 0.03;
	}
	problem->x = nodes;
	problem->l = count;
	problem->y = labels;

	SVMWrapper::createRandomPartitions(problem, 2, partitions);
	TEST_EQUAL(partitions.size(), 2)
	TEST_EQUAL(partitions[0]->l, 2)
	TEST_EQUAL(partitions[1]->l, 2)
END_SECTION

START_SECTION((static void createRandomPartitions(const SVMData &problem, Size number, std::vector< SVMData > &problems)))
	SVMData problem;
	UInt count = 4;
	vector<double> labels;
	vector< vector<pair<Int, double> > > sequences;
	vector<pair<Int, double> > sequence;
	std::vector< SVMData > partitions;
	
	for (Size i = 0; i < count; i++)
	{
		sequence.clear();
		sequence.push_back(make_pair(1, rand()));
		sequence.push_back(make_pair(2, rand()));
		sequence.push_back(make_pair(3, rand()));
		sequence.push_back(make_pair(4, rand()));
		sequences.push_back(sequence);
		labels.push_back(i * 2 / 3 + 0.03);
	}
	problem.sequences = sequences;
	problem.labels = labels;

	SVMWrapper::createRandomPartitions(problem, 2, partitions);
	TEST_EQUAL(partitions.size(), 2)
	TEST_EQUAL(partitions[0].sequences.size(), 2)
	TEST_EQUAL(partitions[1].sequences.size(), 2)
END_SECTION

START_SECTION((static svm_problem* mergePartitions(const std::vector< svm_problem * > &problems, Size except)))
	 svm_problem* problem = new svm_problem();
	 svm_problem* problem2;
	 UInt count = 10;
	 UInt number_of_partitions = 5;
	 svm_node** nodes = new svm_node*[count];
	 double* labels = new double[count];
	 std::vector<svm_problem*> partitions;


	for (Size i = 0; i < count; i++)
	{
		nodes[i] = new svm_node[count];
		nodes[i][count - 1].index = -1;
		labels[i] = ((double) i * 2) / 3 + 0.03;
		for (Size j = 0; j < count; j++)
		{
			nodes[i][j].value = ((double) i * 2) / 3;
		}
	}
	problem->x = nodes;
	problem->l = count;
	problem->y = labels;

	SVMWrapper::createRandomPartitions(problem, number_of_partitions, partitions);
	problem2 = SVMWrapper::mergePartitions(partitions, 4);
	UInt problem2_size = (count / number_of_partitions) * (number_of_partitions - 1);
	UInt partition_size = count / number_of_partitions;
	TEST_EQUAL((UInt) problem2->l, problem2_size)
	for (Size i = 0; i < problem2_size; i++)
	{
		UInt j = 0;
		while(problem->x[i][j].index != -1 && problem2->x[i][j].index != -1)
		{
			TEST_REAL_SIMILAR(partitions[i / partition_size]->x[i % partition_size][j].value, problem2->x[i][j].value)
			++j;
		}
		TEST_REAL_SIMILAR(partitions[i / partition_size]->y[i % partition_size], problem2->y[i])
	}
END_SECTION

START_SECTION((static void mergePartitions(const std::vector< SVMData > &problems, Size except, SVMData &merged_problem)))
	SVMData problem;
	SVMData problem2;
	UInt count = 10;
	UInt number_of_partitions = 5;
	vector<double> labels;
	SVMData merged_problem;
	vector<pair<Int, double> > temp_vector;
	vector<vector<pair<Int, double> > > vectors;
	std::vector< SVMData > partitions;
		
	for (Int i = 0; (UInt) i < count; i++)
	{
		temp_vector.clear();
		labels.push_back(((double) i * 2) / 3 + 0.03);
		for (Int j = 0; (UInt) j < count; j++)
		{
			temp_vector.push_back(make_pair(i * 2, ((double) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	problem.sequences = vectors;
	problem.labels = labels;

	SVMWrapper::createRandomPartitions(problem, number_of_partitions, partitions);
	SVMWrapper::mergePartitions(partitions, 4, problem2);
	UInt problem2_size = (count / number_of_partitions) * (number_of_partitions - 1);
	UInt partition_size = count / number_of_partitions;
	TEST_EQUAL((UInt) problem2.sequences.size(), problem2_size)
	for (UInt i = 0; i < problem2_size; i++)
	{
		UInt j = 0;
		while(j < partitions[i / partition_size].sequences[i % partition_size].size() && j < problem2.sequences[i].size())
		{
			TEST_REAL_SIMILAR(partitions[i / partition_size].sequences[i % partition_size][j].second, problem2.sequences[i][j].second)
			++j;
		}
		TEST_REAL_SIMILAR(partitions[i / partition_size].labels[i % partition_size], problem2.labels[i])
	}
END_SECTION

START_SECTION((static void calculateGaussTable(Size border_length, double sigma, std::vector< double > &gauss_table)))
  UInt border_length = 5;
  double sigma = 2;
  double sigma_square = sigma * sigma;
  vector<double> gauss_table;
  svm.calculateGaussTable(border_length, sigma, gauss_table);

  TEST_EQUAL(gauss_table.size(), 5)
  TEST_EQUAL(gauss_table[0], 1)
  TEST_REAL_SIMILAR(gauss_table[1], exp((-1 / (4.0 * sigma_square)) * 1))
  TEST_REAL_SIMILAR(gauss_table[2], exp((-1 / (4.0 * sigma_square)) * 4))
  TEST_REAL_SIMILAR(gauss_table[3], exp((-1 / (4.0 * sigma_square)) * 9))
  TEST_REAL_SIMILAR(gauss_table[4], exp((-1 / (4.0 * sigma_square)) * 16))
END_SECTION

START_SECTION((double performCrossValidation(svm_problem *problem_ul, const SVMData &problem_l, const bool is_labeled, const std::map< SVM_parameter_type, double > &start_values_map, const std::map< SVM_parameter_type, double > &step_sizes_map, const std::map< SVM_parameter_type, double > &end_values_map, Size number_of_partitions, Size number_of_runs, std::map< SVM_parameter_type, double > &best_parameters, bool additive_step_sizes=true, bool output=false, String performances_file_name="performances.txt", bool mcc_as_performance_measure=false) ))
 
{
  map<SVMWrapper::SVM_parameter_type, double> start_values;
	map<SVMWrapper::SVM_parameter_type, double> step_sizes;
	map<SVMWrapper::SVM_parameter_type, double> end_values;
	LibSVMEncoder encoder;
	vector< vector< pair<Int, double> > > vectors;
	vector< pair<Int, double> > temp_vector;
	vector<svm_node*> encoded_vectors;
	UInt count = 8;
	vector<double> labels;
	svm_problem* problem;
	map<SVMWrapper::SVM_parameter_type, double> parameters;
	double cv_quality;

	for (UInt j = 0; j < count; j++)
	{
		temp_vector.clear();
		for (UInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((double) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoder.encodeLibSVMVectors(vectors, encoded_vectors);
	for (Size i = 0; i < count; i++)
	{
		labels.push_back(((double) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);

	start_values.insert(make_pair(SVMWrapper::C, 1));
	step_sizes.insert(make_pair(SVMWrapper::C, 100));
	end_values.insert(make_pair(SVMWrapper::C, 1000));

	start_values.insert(make_pair(SVMWrapper::NU, 0.4));
	step_sizes.insert(make_pair(SVMWrapper::NU, 0.1));
	end_values.insert(make_pair(SVMWrapper::NU, 0.6));

	start_values.insert(make_pair(SVMWrapper::DEGREE, 1));
	step_sizes.insert(make_pair(SVMWrapper::DEGREE, 1));
	end_values.insert(make_pair(SVMWrapper::DEGREE, 3));

  SVMData problem_2;
	cv_quality = svm.performCrossValidation(problem, problem_2, false, start_values, step_sizes, end_values, 2, 1, parameters, true, false);
	TEST_NOT_EQUAL(parameters.size(), 0)
  TEST_REAL_SIMILAR(cv_quality, 1)
}
  // CV, method 2
    
  map<SVMWrapper::SVM_parameter_type, double> start_values;
	map<SVMWrapper::SVM_parameter_type, double> step_sizes;
	map<SVMWrapper::SVM_parameter_type, double> end_values;
	LibSVMEncoder encoder;
	UInt count = 8;
	map<SVMWrapper::SVM_parameter_type, double> parameters;
	double cv_quality;
	SVMWrapper svm2;
	SVMData problem;
	vector<double> labels;
	vector< vector<pair<Int, double> > > sequences;
	vector<pair<Int, double> > sequence;
	
	svm2.setParameter(SVMWrapper::KERNEL_TYPE, SVMWrapper::OLIGO);
	svm2.setParameter(SVMWrapper::BORDER_LENGTH, 2);
	svm2.setParameter(SVMWrapper::C, 1);
	svm2.setParameter(SVMWrapper::SIGMA, 1);
	svm2.setParameter(SVMWrapper::SVM_TYPE, NU_SVR);

	for (Size i = 0; i < count; i++)
	{
		sequence.clear();
		sequence.push_back(make_pair(1, rand()));
		sequence.push_back(make_pair(2, rand()));
		sequence.push_back(make_pair(3, rand()));
		sequence.push_back(make_pair(4, rand()));
		sequences.push_back(sequence);
		labels.push_back(i * 2 / 3 + 0.03);
	}
	problem.sequences = sequences;
	problem.labels = labels;

	start_values.insert(make_pair(SVMWrapper::C, 1));
	step_sizes.insert(make_pair(SVMWrapper::C, 100));
	end_values.insert(make_pair(SVMWrapper::C, 1000));

	start_values.insert(make_pair(SVMWrapper::NU, 0.4));
	step_sizes.insert(make_pair(SVMWrapper::NU, 0.1));
	end_values.insert(make_pair(SVMWrapper::NU, 0.6));

	start_values.insert(make_pair(SVMWrapper::DEGREE, 1));
	step_sizes.insert(make_pair(SVMWrapper::DEGREE, 1));
	end_values.insert(make_pair(SVMWrapper::DEGREE, 3));

  cv_quality = svm2.performCrossValidation(nullptr, problem, true, start_values, step_sizes, end_values, 2, 1, parameters, true, false);

  TEST_NOT_EQUAL(parameters.size(), 0)
  // cv_quality is nan
  TEST_EQUAL(cv_quality != cv_quality, true)
END_SECTION

START_SECTION((void predict(struct svm_problem *problem, std::vector< double > &predicted_labels)))
 	LibSVMEncoder encoder;
	vector< vector< pair<Int, double> > > vectors;
	vector< pair<Int, double> > temp_vector;
	vector<svm_node*> encoded_vectors;
	UInt count = 8;
	vector<double> labels;
	vector<double> predicted_labels;
	svm_problem* problem;

	for (UInt j = 0; j < count; j++)
	{
		temp_vector.clear();
		for (UInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((double) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoder.encodeLibSVMVectors(vectors, encoded_vectors);
	for (Size i = 0; i < count; i++)
	{
		labels.push_back(((double) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.train(problem);
	svm.predict(problem, predicted_labels);
	TEST_NOT_EQUAL(predicted_labels.size(), 0)
END_SECTION

START_SECTION((void predict(const SVMData &problem, std::vector< double > &results)))
	SVMWrapper svm2;
 	LibSVMEncoder encoder;
	vector< vector< pair<Int, double> > > sequences;
	vector< pair<Int, double> > sequence;
	UInt count = 8;
	vector<double> labels;
	vector<double> predicted_labels;
	SVMData problem;

	svm2.setParameter(SVMWrapper::KERNEL_TYPE, SVMWrapper::OLIGO);
	svm2.setParameter(SVMWrapper::BORDER_LENGTH, 2);
	svm2.setParameter(SVMWrapper::C, 1);
	svm2.setParameter(SVMWrapper::SIGMA, 1);
	svm2.setParameter(SVMWrapper::SVM_TYPE, NU_SVR);

	for (Size i = 0; i < count; i++)
	{
		sequence.clear();
		sequence.push_back(make_pair(1, rand()));
		sequence.push_back(make_pair(2, rand()));
		sequence.push_back(make_pair(3, rand()));
		sequence.push_back(make_pair(4, rand()));
		sequences.push_back(sequence);
		labels.push_back(i * 2 / 3 + 0.03);
	}

	problem.sequences = sequences;
	problem.labels = labels;
	svm2.train(problem);
	svm2.predict(problem, predicted_labels);
	TEST_NOT_EQUAL(predicted_labels.size(), 0)
END_SECTION

START_SECTION((svm_problem* computeKernelMatrix(svm_problem* problem1, svm_problem* problem2)))
	vector<String> sequences;
	String allowed_characters = "ACNGT";
	String output;
	Int border_length = 5;
	double sigma = 2;
	vector< pair<Int, double> > encoded_sequence;
  vector<double> labels;
  struct svm_problem* data;
  struct svm_problem* kernel_matrix;
  LibSVMEncoder encoder;

	svm.setParameter(SVMWrapper::BORDER_LENGTH, border_length);
	svm.setParameter(SVMWrapper::SIGMA, sigma);
	svm.setParameter(SVMWrapper::KERNEL_TYPE, SVMWrapper::OLIGO);
  labels.push_back(1);
  labels.push_back(2);
  sequences.push_back("ACNNGTATCA");
  sequences.push_back("AACNNGTACCA");
	data = encoder.encodeLibSVMProblemWithOligoBorderVectors(sequences, labels, 1, allowed_characters, border_length);
	kernel_matrix = svm.computeKernelMatrix(data, data);
	svm.train(data);

	TOLERANCE_ABSOLUTE(0.0001)
	TEST_REAL_SIMILAR(kernel_matrix->x[0][0].value, 1)
	TEST_REAL_SIMILAR(kernel_matrix->x[0][1].value, 19.7156)
	TEST_REAL_SIMILAR(kernel_matrix->x[0][2].value, 21.1308)
	TEST_REAL_SIMILAR(kernel_matrix->x[1][0].value, 2)
	TEST_REAL_SIMILAR(kernel_matrix->x[1][1].value, 21.1308)
	TEST_REAL_SIMILAR(kernel_matrix->x[1][2].value, 27.2309)
	TEST_EQUAL(kernel_matrix->x[0][0].index, 0)
	TEST_EQUAL(kernel_matrix->x[0][1].index, 1)
	TEST_EQUAL(kernel_matrix->x[0][2].index, 2)
	TEST_EQUAL(kernel_matrix->x[1][0].index, 0)
	TEST_EQUAL(kernel_matrix->x[1][1].index, 1)
	TEST_EQUAL(kernel_matrix->x[1][2].index, 2)
	TEST_EQUAL(kernel_matrix->y[0], 1)
	TEST_EQUAL(kernel_matrix->y[1], 2)

END_SECTION

START_SECTION((svm_problem* computeKernelMatrix(const SVMData &problem1, const SVMData &problem2)))
	vector<AASequence> sequences;
	String allowed_characters = "ACNGT";
	String output;
	Int border_length = 5;
	double sigma = 2;
  vector<double> labels;
  struct svm_problem* kernel_matrix;
  LibSVMEncoder encoder;
	vector< vector< std::pair< Int, double > > > data;
	SVMData svm_data;

	svm.setParameter(SVMWrapper::BORDER_LENGTH, border_length);
	svm.setParameter(SVMWrapper::SIGMA, sigma);
	svm.setParameter(SVMWrapper::KERNEL_TYPE, SVMWrapper::OLIGO);
  labels.push_back(1);
  labels.push_back(2);
  sequences.push_back(AASequence::fromString("ACNNGTATCA"));
  sequences.push_back(AASequence::fromString("AACNNGTACCA"));
	encoder.encodeProblemWithOligoBorderVectors(sequences, 1, allowed_characters, border_length, data);
	svm_data.sequences = data;
	svm_data.labels = labels;
	
	kernel_matrix = svm.computeKernelMatrix(svm_data, svm_data);
			
	TOLERANCE_ABSOLUTE(0.0001)
	TEST_REAL_SIMILAR(kernel_matrix->x[0][0].value, 1)
	TEST_REAL_SIMILAR(kernel_matrix->x[0][1].value, 19.7156)
	TEST_REAL_SIMILAR(kernel_matrix->x[0][2].value, 21.1308)
	TEST_REAL_SIMILAR(kernel_matrix->x[1][0].value, 2)
	TEST_REAL_SIMILAR(kernel_matrix->x[1][1].value, 21.1308)
	TEST_REAL_SIMILAR(kernel_matrix->x[1][2].value, 27.2309)
	TEST_EQUAL(kernel_matrix->x[0][0].index, 0)
	TEST_EQUAL(kernel_matrix->x[0][1].index, 1)
	TEST_EQUAL(kernel_matrix->x[0][2].index, 2)
	TEST_EQUAL(kernel_matrix->x[1][0].index, 0)
	TEST_EQUAL(kernel_matrix->x[1][1].index, 1)
	TEST_EQUAL(kernel_matrix->x[1][2].index, 2)
	TEST_EQUAL(kernel_matrix->y[0], 1)
	TEST_EQUAL(kernel_matrix->y[1], 2)

END_SECTION

START_SECTION((static double kernelOligo(const svm_node *x, const svm_node *y, const std::vector< double > &gauss_table, double sigma_square=0, Size max_distance=50)))
  vector<double> labels;
	String sequence = "ACNNGTATCA";
	String allowed_characters = "ACNGT";
	String output;
	Int border_length = 5;
	svm_problem* data;
	double result = 0;
  double sigma = 2;
  vector<double> gauss_table;
	vector<String> sequences;
  svm.calculateGaussTable(border_length, sigma, gauss_table);
	LibSVMEncoder encoder;
	svm.setParameter(SVMWrapper::BORDER_LENGTH, border_length);
	svm.setParameter(SVMWrapper::SIGMA, sigma);
	svm.setParameter(SVMWrapper::KERNEL_TYPE, SVMWrapper::OLIGO);

  labels.push_back(1);
  labels.push_back(2);
  sequences.push_back("ACNNGTATCA");
  sequences.push_back("AACNNGTACCA");
	data = encoder.encodeLibSVMProblemWithOligoBorderVectors(sequences, labels, 1, allowed_characters, border_length);
	result = SVMWrapper::kernelOligo(data->x[0], data->x[1], gauss_table);
	TOLERANCE_ABSOLUTE(0.0001)
	TEST_REAL_SIMILAR(result, 21.1308)
	delete data;
END_SECTION

START_SECTION((static double kernelOligo(const std::vector< std::pair< int, double > > &x, const std::vector< std::pair< int, double > > &y, const std::vector< double > &gauss_table, int max_distance=-1)))
  vector<double> labels;
	String sequence = "ACNNGTATCA";
	String allowed_characters = "ACNGT";
	String output;
	Int border_length = 5;
	vector< vector< std::pair< Int, double > > > data;
	double result = 0;
  double sigma = 2;
  vector<double> gauss_table;
	vector<AASequence> sequences;
  svm.calculateGaussTable(border_length, sigma, gauss_table);
	LibSVMEncoder encoder;
	svm.setParameter(SVMWrapper::BORDER_LENGTH, border_length);
	svm.setParameter(SVMWrapper::SIGMA, sigma);
	svm.setParameter(SVMWrapper::KERNEL_TYPE, SVMWrapper::OLIGO);

  labels.push_back(1);
  labels.push_back(2);
  sequences.push_back(AASequence::fromString("ACNNGTATCA"));
  sequences.push_back(AASequence::fromString("AACNNGTACCA"));
	encoder.encodeProblemWithOligoBorderVectors(sequences, 1, allowed_characters, border_length, data);
	result = SVMWrapper::kernelOligo(data[0], data[1], gauss_table);
	TOLERANCE_ABSOLUTE(0.0001)
	TEST_REAL_SIMILAR(result, 21.1308)
END_SECTION

START_SECTION((void getDecisionValues(svm_problem* data, std::vector<double>& decision_values)))
 	LibSVMEncoder encoder;
	vector< vector< pair<Int, double> > > vectors;
	vector< pair<Int, double> > temp_vector;
	vector<svm_node*> encoded_vectors;
	UInt count = 8;
	vector<double> labels;
	vector<double> predicted_labels;
	svm_problem* problem;
	vector<double> decision_values;

	svm.setParameter(SVMWrapper::SVM_TYPE, NU_SVR);
	svm.setParameter(SVMWrapper::KERNEL_TYPE, POLY);
	svm.setParameter(SVMWrapper::DEGREE, 2);
	for (UInt j = 0; j < count; j++)
	{
		temp_vector.clear();
		for (UInt i = 1; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((double) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoder.encodeLibSVMVectors(vectors, encoded_vectors);
	for (Size i = 0; i < count; i++)
	{
		labels.push_back(((double) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.train(problem);
	svm.predict(problem, predicted_labels);
	TEST_NOT_EQUAL(predicted_labels.size(), 0)
	svm.getDecisionValues(problem, decision_values);
	TEST_EQUAL(predicted_labels == decision_values, true)

	svm.setParameter(SVMWrapper::SVM_TYPE, C_SVC);
	labels.clear();
	labels.resize(4, 1);
	labels.resize(8, -1);
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.train(problem);
	svm.predict(problem, predicted_labels);
	TEST_NOT_EQUAL(predicted_labels.size(), 0)
	svm.getDecisionValues(problem, decision_values);
	TEST_EQUAL(predicted_labels.size() == decision_values.size(), true)
	for (Size i = 0; i < predicted_labels.size(); ++i)
	{
		TEST_EQUAL((predicted_labels[i] < 0 && decision_values[i] < 0)
							|| (predicted_labels[i] > 0 && decision_values[i] > 0), true)
	}
	labels.clear();
	labels.resize(4, -1);
	labels.resize(8, 1);
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.train(problem);
	svm.predict(problem, predicted_labels);
	TEST_NOT_EQUAL(predicted_labels.size(), 0)
	svm.getDecisionValues(problem, decision_values);
	TEST_EQUAL(predicted_labels.size() == decision_values.size(), true)
	for (Size i = 0; i < predicted_labels.size(); ++i)
	{
		TEST_EQUAL((predicted_labels[i] < 0 && decision_values[i] < 0)
							|| (predicted_labels[i] > 0 && decision_values[i] > 0), true)
	}

END_SECTION

START_SECTION((void scaleData(svm_problem* data, Int max_scale_value = -1)))
 	LibSVMEncoder encoder;
	vector< vector< pair<Int, double> > > vectors;
	vector< pair<Int, double> > temp_vector;
	vector<svm_node*> encoded_vectors;
	UInt count = 8;
	vector<double> labels;
	svm_problem* problem;
	vector<double> decision_values;

	svm.setParameter(SVMWrapper::SVM_TYPE, NU_SVR);
	svm.setParameter(SVMWrapper::KERNEL_TYPE, POLY);
	svm.setParameter(SVMWrapper::DEGREE, 2);
	for (UInt j = 0; j < count; j++)
	{
		temp_vector.clear();
		for (UInt i = 1; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i, ((double) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoder.encodeLibSVMVectors(vectors, encoded_vectors);
	for (Size i = 0; i < count; i++)
	{
		labels.push_back(((double) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.scaleData(problem, 2);

	TEST_REAL_SIMILAR(problem->x[0][0].value, 0)
	TEST_REAL_SIMILAR(problem->x[0][1].value, 0)
	TEST_REAL_SIMILAR(problem->x[0][2].value, 0)
	TEST_REAL_SIMILAR(problem->x[0][3].value, 0)
	TEST_REAL_SIMILAR(problem->x[0][4].value, 0)
	TEST_REAL_SIMILAR(problem->x[1][0].value, 0.2857)
	TEST_REAL_SIMILAR(problem->x[1][1].value, 0.2857)
	TEST_REAL_SIMILAR(problem->x[1][2].value, 0.2857)
	TEST_REAL_SIMILAR(problem->x[1][3].value, 0.2857)
	TEST_REAL_SIMILAR(problem->x[1][4].value, 0.2857)
	TEST_REAL_SIMILAR(problem->x[2][0].value, 0.5714)
	TEST_REAL_SIMILAR(problem->x[2][1].value, 0.5714)
	TEST_REAL_SIMILAR(problem->x[2][2].value, 0.5714)
	TEST_REAL_SIMILAR(problem->x[2][3].value, 0.5714)
	TEST_REAL_SIMILAR(problem->x[2][4].value, 0.5714)
	TEST_REAL_SIMILAR(problem->x[3][0].value, 0.8571)
	TEST_REAL_SIMILAR(problem->x[3][1].value, 0.8571)
	TEST_REAL_SIMILAR(problem->x[3][2].value, 0.8571)
	TEST_REAL_SIMILAR(problem->x[3][3].value, 0.8571)
	TEST_REAL_SIMILAR(problem->x[3][4].value, 0.8571)
	TEST_REAL_SIMILAR(problem->x[4][0].value, 1.1429)
	TEST_REAL_SIMILAR(problem->x[4][1].value, 1.1429)
	TEST_REAL_SIMILAR(problem->x[4][2].value, 1.1429)
	TEST_REAL_SIMILAR(problem->x[4][3].value, 1.1429)
	TEST_REAL_SIMILAR(problem->x[4][4].value, 1.1429)
	TEST_REAL_SIMILAR(problem->x[5][0].value, 1.4286)
	TEST_REAL_SIMILAR(problem->x[5][1].value, 1.4286)
	TEST_REAL_SIMILAR(problem->x[5][2].value, 1.4286)
	TEST_REAL_SIMILAR(problem->x[5][3].value, 1.4286)
	TEST_REAL_SIMILAR(problem->x[5][4].value, 1.4286)
	TEST_REAL_SIMILAR(problem->x[6][0].value, 1.7143)
	TEST_REAL_SIMILAR(problem->x[6][1].value, 1.7143)
	TEST_REAL_SIMILAR(problem->x[6][2].value, 1.7143)
	TEST_REAL_SIMILAR(problem->x[6][3].value, 1.7143)
	TEST_REAL_SIMILAR(problem->x[6][4].value, 1.7143)
	TEST_REAL_SIMILAR(problem->x[7][0].value, 2)
	TEST_REAL_SIMILAR(problem->x[7][1].value, 2)
	TEST_REAL_SIMILAR(problem->x[7][2].value, 2)
	TEST_REAL_SIMILAR(problem->x[7][3].value, 2)
	TEST_REAL_SIMILAR(problem->x[7][4].value, 2)

	svm.scaleData(problem);

	TEST_REAL_SIMILAR(problem->x[0][0].value, -1)
	TEST_REAL_SIMILAR(problem->x[0][1].value, -1)
	TEST_REAL_SIMILAR(problem->x[0][2].value, -1)
	TEST_REAL_SIMILAR(problem->x[0][3].value, -1)
	TEST_REAL_SIMILAR(problem->x[0][4].value, -1)
	TEST_REAL_SIMILAR(problem->x[1][0].value, -0.7143)
	TEST_REAL_SIMILAR(problem->x[1][1].value, -0.7143)
	TEST_REAL_SIMILAR(problem->x[1][2].value, -0.7143)
	TEST_REAL_SIMILAR(problem->x[1][3].value, -0.7143)
	TEST_REAL_SIMILAR(problem->x[1][4].value, -0.7143)
	TEST_REAL_SIMILAR(problem->x[2][0].value, -0.4286)
	TEST_REAL_SIMILAR(problem->x[2][1].value, -0.4286)
	TEST_REAL_SIMILAR(problem->x[2][2].value, -0.4286)
	TEST_REAL_SIMILAR(problem->x[2][3].value, -0.4286)
	TEST_REAL_SIMILAR(problem->x[2][4].value, -0.4286)
	TEST_REAL_SIMILAR(problem->x[3][0].value, -0.1429)
	TEST_REAL_SIMILAR(problem->x[3][1].value, -0.1429)
	TEST_REAL_SIMILAR(problem->x[3][2].value, -0.1429)
	TEST_REAL_SIMILAR(problem->x[3][3].value, -0.1429)
	TEST_REAL_SIMILAR(problem->x[3][4].value, -0.1429)
	TEST_REAL_SIMILAR(problem->x[4][0].value, 0.1429)
	TEST_REAL_SIMILAR(problem->x[4][1].value, 0.1429)
	TEST_REAL_SIMILAR(problem->x[4][2].value, 0.1429)
	TEST_REAL_SIMILAR(problem->x[4][3].value, 0.1429)
	TEST_REAL_SIMILAR(problem->x[4][4].value, 0.1429)
	TEST_REAL_SIMILAR(problem->x[5][0].value, 0.4286)
	TEST_REAL_SIMILAR(problem->x[5][1].value, 0.4286)
	TEST_REAL_SIMILAR(problem->x[5][2].value, 0.4286)
	TEST_REAL_SIMILAR(problem->x[5][3].value, 0.4286)
	TEST_REAL_SIMILAR(problem->x[5][4].value, 0.4286)
	TEST_REAL_SIMILAR(problem->x[6][0].value, 0.7143)
	TEST_REAL_SIMILAR(problem->x[6][1].value, 0.7143)
	TEST_REAL_SIMILAR(problem->x[6][2].value, 0.7143)
	TEST_REAL_SIMILAR(problem->x[6][3].value, 0.7143)
	TEST_REAL_SIMILAR(problem->x[6][4].value, 0.7143)
	TEST_REAL_SIMILAR(problem->x[7][0].value, 1)
	TEST_REAL_SIMILAR(problem->x[7][1].value, 1)
	TEST_REAL_SIMILAR(problem->x[7][2].value, 1)
	TEST_REAL_SIMILAR(problem->x[7][3].value, 1)
	TEST_REAL_SIMILAR(problem->x[7][4].value, 1)

END_SECTION

START_SECTION((void getSignificanceBorders(svm_problem *data, std::pair< double, double > &borders, double confidence=0.95, Size number_of_runs=5, Size number_of_partitions=5, double step_size=0.01, Size max_iterations=1000000)))
	NOT_TESTABLE
END_SECTION

START_SECTION((void getSignificanceBorders(const SVMData &data, std::pair< double, double > &sigmas, double confidence=0.95, Size number_of_runs=5, Size number_of_partitions=5, double step_size=0.01, Size max_iterations=1000000)))
	NOT_TESTABLE
END_SECTION

START_SECTION((double getPValue(double sigma1, double sigma2, std::pair<double, double> point)))

	pair<double, double> point;

	point.first = 0.447934;
	point.second = 0.404208;

	TEST_REAL_SIMILAR(svm.getPValue(0.18, 1.06, point), 0.327505)
END_SECTION

START_SECTION((void setTrainingSample(svm_problem* training_sample)))
	NOT_TESTABLE
END_SECTION

START_SECTION((void setTrainingSample(SVMData &training_sample)))
	NOT_TESTABLE
END_SECTION

START_SECTION((void setParameter(SVM_parameter_type type, double value)))
 	svm.setParameter(SVMWrapper::C, 1.0043);
	svm.setParameter(SVMWrapper::NU, 0.0523);
	svm.setParameter(SVMWrapper::P, 1.2319);

	TEST_REAL_SIMILAR(svm.getDoubleParameter(SVMWrapper::C), 1.0043)
	TEST_REAL_SIMILAR(svm.getDoubleParameter(SVMWrapper::NU), 0.0523)
	TEST_REAL_SIMILAR(svm.getDoubleParameter(SVMWrapper::P), 1.2319)
END_SECTION

START_SECTION((void setParameter(SVM_parameter_type type, Int value)))
	svm.setParameter(SVMWrapper::SVM_TYPE, EPSILON_SVR);
	svm.setParameter(SVMWrapper::KERNEL_TYPE, LINEAR);
	svm.setParameter(SVMWrapper::DEGREE, 2);
	svm.setParameter(SVMWrapper::C, 23);
	svm.setParameter(SVMWrapper::PROBABILITY, 1);

	TEST_EQUAL(svm.getIntParameter(SVMWrapper::SVM_TYPE)==EPSILON_SVR,true);
	TEST_EQUAL(svm.getIntParameter(SVMWrapper::KERNEL_TYPE)==LINEAR,true);
	TEST_EQUAL(svm.getIntParameter(SVMWrapper::DEGREE)==2,true);
	TEST_EQUAL((int) svm.getDoubleParameter(SVMWrapper::C), 23);
	TEST_EQUAL(svm.getIntParameter(SVMWrapper::PROBABILITY), 1)
END_SECTION

START_SECTION((virtual ~SVMWrapper()))
	delete ptr;
END_SECTION

START_SECTION((void loadModel(std::string modelFilename)))
	LibSVMEncoder encoder;
	svm.setParameter(SVMWrapper::KERNEL_TYPE, POLY);
	vector< vector< pair<Int, double> > > vectors;
	vector< pair<Int, double> > temp_vector;
	vector<svm_node*> encoded_vectors;
	UInt count = 8;
	vector<double> labels;
	vector<double> predicted_labels1;
	vector<double> predicted_labels2;
	svm_problem* problem;
	SVMWrapper svm2;

	for (UInt j = 0; j < count; j++)
	{
		temp_vector.clear();
		for (UInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((double) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoder.encodeLibSVMVectors(vectors, encoded_vectors);
	for (Size i = 0; i < count; i++)
	{
		labels.push_back(((double) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.train(problem);
	svm.predict(problem, predicted_labels1);

	String filename = "svm.model";
	NEW_TMP_FILE(filename)
	svm.saveModel(filename);
	svm2.loadModel(filename);
	svm2.predict(problem, predicted_labels2);
	TEST_NOT_EQUAL(predicted_labels1.size(), 0)
	TEST_EQUAL(predicted_labels1.size(), predicted_labels2.size())
	for (Size i = 0; i < predicted_labels1.size(); i++)
	{
		TEST_REAL_SIMILAR(predicted_labels1[i], predicted_labels2[i])
	}
END_SECTION

START_SECTION((void saveModel(std::string modelFilename) const))
	LibSVMEncoder encoder;
	svm.setParameter(SVMWrapper::KERNEL_TYPE, POLY);
	vector< vector< pair<Int, double> > > vectors;
	vector< pair<Int, double> > temp_vector;
	vector<svm_node*> encoded_vectors;
	UInt count = 8;
	vector<double> labels;
	vector<double> predicted_labels1;
	vector<double> predicted_labels2;
	svm_problem* problem;
	SVMWrapper svm2;

	for (UInt j = 0; j < count; j++)
	{
		temp_vector.clear();
		for (UInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((double) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoder.encodeLibSVMVectors(vectors, encoded_vectors);
	for (Size i = 0; i < count; i++)
	{
		labels.push_back(((double) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.train(problem);

	String filename = "svm.model";
	NEW_TMP_FILE(filename)
	svm.saveModel(filename);
	svm2.loadModel(filename);
	svm.predict(problem, predicted_labels1);
	svm2.predict(problem, predicted_labels2);
	TEST_NOT_EQUAL(predicted_labels1.size(), 0)
	TEST_NOT_EQUAL(predicted_labels2.size(), 0)
	TEST_EQUAL(predicted_labels1.size(), predicted_labels2.size())

	for (Size i = 0; i < predicted_labels1.size(); i++)
	{
		TEST_REAL_SIMILAR(predicted_labels1[i], predicted_labels2[i])
	}
END_SECTION

START_SECTION((void predict(const std::vector< svm_node * > &vectors, std::vector< double > &predicted_rts)))
	LibSVMEncoder encoder;
	vector< vector< pair<Int, double> > > vectors;
	vector< pair<Int, double> > temp_vector;
	vector<svm_node*> encoded_vectors;
	UInt count = 8;
	vector<double> labels;
	vector<double> predicted_labels;
	svm_problem* problem;

	for (UInt j = 0; j < count; j++)
	{
		for (UInt i = 0; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((double) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoder.encodeLibSVMVectors(vectors, encoded_vectors);
	for (Size i = 0; i < count; i++)
	{
		labels.push_back(((double) i * 2) / 3 + 0.03);
	}
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.train(problem);
	svm.predict(encoded_vectors, predicted_labels);
	TEST_NOT_EQUAL(predicted_labels.size(), 0)

END_SECTION

START_SECTION((void setWeights(const std::vector< Int > &weight_labels, const std::vector< double > &weights)))
	NOT_TESTABLE
END_SECTION

START_SECTION((void getSVCProbabilities(struct svm_problem *problem, std::vector< double > &probabilities, std::vector< double > &prediction_labels)))
 	LibSVMEncoder encoder;
	vector< vector< pair<Int, double> > > vectors;
	vector< pair<Int, double> > temp_vector;
	vector<svm_node*> encoded_vectors;
	UInt count = 8;
	vector<double> labels;
	vector<double> predicted_labels;
	svm_problem* problem;
	vector<double> probabilities;

	svm.setParameter(SVMWrapper::SVM_TYPE, C_SVC);
	svm.setParameter(SVMWrapper::KERNEL_TYPE, POLY);
	svm.setParameter(SVMWrapper::DEGREE, 2);
	for (UInt j = 0; j < count; j++)
	{
		temp_vector.clear();
		for (UInt i = 1; i < 6; i++)
		{
			temp_vector.push_back(make_pair(i * 2, ((double) i) * j * 0.3));
		}
		vectors.push_back(temp_vector);
	}
	encoder.encodeLibSVMVectors(vectors, encoded_vectors);

	labels.clear();
	labels.resize(count / 2, 1);
	labels.resize(count, -1);
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.train(problem);
	svm.predict(problem, predicted_labels);
	TEST_NOT_EQUAL(predicted_labels.size(), 0)
	svm.getSVCProbabilities(problem, probabilities, predicted_labels);
	TEST_EQUAL(predicted_labels.size() == probabilities.size(), true)
	for (Size i = 0; i < predicted_labels.size(); ++i)
	{
		TEST_EQUAL((predicted_labels[i] < 0 && probabilities[i] < 0.5)
							|| (predicted_labels[i] > 0 && probabilities[i] >= 0.5), true)
	}
	labels.clear();
	// Start with -1 as "first" label
	labels.resize(count / 2, -1);
	labels.resize(count, 1);
	problem = encoder.encodeLibSVMProblem(encoded_vectors, labels);
	svm.train(problem);
	svm.predict(problem, predicted_labels);
	TEST_NOT_EQUAL(predicted_labels.size(), 0)
	svm.getSVCProbabilities(problem, probabilities, predicted_labels);
	TEST_EQUAL(predicted_labels.size() == probabilities.size(), true)
	for (Size i = 0; i < predicted_labels.size(); ++i)
	{
		// At probability 0.5, LibSVM will assign the first encountered label in the training data
		// in this case "-1"
		TEST_EQUAL((predicted_labels[i] < 0 && probabilities[i] <= 0.5)
							|| (predicted_labels[i] > 0 && probabilities[i] > 0.5), true)
	}

END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
