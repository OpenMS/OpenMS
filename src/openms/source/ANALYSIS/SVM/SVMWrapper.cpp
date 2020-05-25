// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/CONCEPT/LogStream.h>


#include <fstream>

#include <boost/math/distributions/normal.hpp>

using namespace std;
using boost::math::cdf;

namespace OpenMS
{

  SVMData::SVMData() :
    sequences(std::vector<std::vector<std::pair<Int, double> > >()),
    labels(std::vector<double>())
  {
  }

  SVMData::SVMData(std::vector<std::vector<std::pair<Int, double> > >& seqs, std::vector<double>& lbls) :
    sequences(seqs),
    labels(lbls)
  {
  }

  bool SVMData::operator==(const SVMData& rhs) const
  {
    return sequences == rhs.sequences
           && labels == rhs.labels;
  }

  bool SVMData::store(const String& filename) const
  {
    std::ofstream output_file(filename.c_str());

    // checking if file is writable
    if (!File::writable(filename) || sequences.size() != labels.size())
    {
      return false;
    }

    // writing feature vectors
    for (Size i = 0; i < sequences.size(); i++)
    {
      output_file << labels[i] << " ";
      for (Size j = 0; j < sequences[i].size(); ++j)
      {
        output_file << sequences[i][j].second << ":" << sequences[i][j].first << " ";
      }
      output_file << std::endl;
    }
    output_file.flush();
    output_file.close();
    std::cout.flush();
    return true;
  }

  bool SVMData::load(const String& filename)
  {
    Size counter = 0;
    std::vector<String> parts;
    std::vector<String> temp_parts;

    if (!File::exists(filename))
    {
      return false;
    }
    if (!File::readable(filename))
    {
      return false;
    }
    if (File::empty(filename))
    {
      return false;
    }

    TextFile text_file(filename.c_str(), true);
    TextFile::ConstIterator it;

    it = text_file.begin();

    sequences.resize(text_file.end() - text_file.begin(), std::vector<std::pair<Int, double> >());
    labels.resize(text_file.end() - text_file.begin(), 0.);
    while (it != text_file.end())
    {
      it->split(' ', parts);
      labels[counter] = parts[0].trim().toFloat();
      sequences[counter].resize(parts.size(), std::pair<Int, double>());
      for (Size j = 1; j < parts.size(); ++j)
      {
        parts[j].split(':', temp_parts);
        if (temp_parts.size() < 2)
        {
          return false;
        }
        sequences[counter][j - 1].second = temp_parts[0].trim().toFloat();
        sequences[counter][j - 1].first = temp_parts[1].trim().toInt();
      }
      ++counter;
      ++it;
    }
    return true;
  }

  SVMWrapper::SVMWrapper() :
    ProgressLogger(),
    param_(nullptr),
    model_(nullptr),
    sigma_(0),
    sigmas_(vector<double>()),
    gauss_table_(),
    kernel_type_(PRECOMPUTED),
    border_length_(0),
    training_set_(nullptr),
    training_problem_(nullptr),
    training_data_(SVMData())
  {
    param_ = (struct svm_parameter*) malloc(sizeof(struct svm_parameter));
    initParameters_();
  }

  SVMWrapper::~SVMWrapper()
  {
    if (param_ != nullptr)
    {
      svm_destroy_param(param_);
      free(param_);
      param_ = nullptr;
    }
    if (model_ != nullptr)
    {
#if OPENMS_LIBSVM_VERSION_MAJOR == 2
      svm_destroy_model(model_);
#else
      svm_free_and_destroy_model(&model_);
#endif
      model_ = nullptr;
    }
  }

  void SVMWrapper::setParameter(SVM_parameter_type type, Int value)
  {

    switch (type)
    {
    case (SVM_TYPE):
      if (value == NU_SVR || value == EPSILON_SVR || value == NU_SVC || value == C_SVC || value == ONE_CLASS)
      {
        param_->svm_type = value;
      }
      break;

    case (KERNEL_TYPE):
      kernel_type_ = value;
      if (value == OLIGO)
      {
        param_->kernel_type = PRECOMPUTED;
      }
      else
      {
        param_->kernel_type = value;
      }
      break;

    case (DEGREE):
      param_->degree = value;
      break;

    case (BORDER_LENGTH):
      border_length_ = value;
      break;

    case (C):
      param_->C = value;
      break;

    case (P):
      param_->p = value;
      break;

    case (NU):
      param_->nu = value;
      break;

    case (GAMMA):
      param_->gamma = value;
      break;

    case (SIGMA):
      sigma_ = value;
      if (border_length_ >= 1)
      {
        SVMWrapper::calculateGaussTable(border_length_, sigma_, gauss_table_);
      }
      break;

    case (PROBABILITY):
      if (value == 1 || value == 0)
      {
        param_->probability = value;
      }
      break;
    }
  }

  Int SVMWrapper::getIntParameter(SVM_parameter_type type)
  {
    switch (type)
    {
    case (KERNEL_TYPE):
      if (param_->kernel_type != PRECOMPUTED)
      {
        return param_->kernel_type;
      }
      else
      {
        return (Int)kernel_type_;
      }

    case (SVM_TYPE):
      return param_->svm_type;

    case (DEGREE):
      return (Int) param_->degree;

    case (BORDER_LENGTH):
      return (Int)border_length_;

    case (PROBABILITY):
      return param_->probability;

    default:
      return -1;
    }
  }

  void SVMWrapper::setParameter(SVM_parameter_type type, double value)
  {
    switch (type)
    {
    case (DEGREE):
      param_->degree = (int)value;
      break;

    case (C):
      param_->C = value;
      break;

    case (P):
      param_->p = value;
      break;

    case (NU):
      param_->nu = value;
      break;

    case (GAMMA):
      param_->gamma = value;
      break;

    case (SIGMA):
      sigma_ = value;
      if (border_length_ >= 1)
      {
        SVMWrapper::calculateGaussTable(border_length_, sigma_, gauss_table_);
      }
      break;

    default:
      break;
    }
  }

  double SVMWrapper::getDoubleParameter(SVM_parameter_type type)
  {
    switch (type)
    {
    case (C):
      return param_->C;

    case (P):
      return param_->p;

    case (NU):
      return param_->nu;

    case (GAMMA):
      return param_->gamma;

    case (SIGMA):
      return sigma_;

    default:
      return -1;
    }
  }

  void SVMWrapper::setTrainingSample(svm_problem* training_sample)
  {
    training_set_ = training_sample;
  }

  void SVMWrapper::setTrainingSample(SVMData& training_sample)
  {
    training_data_ = training_sample;
  }

  Int SVMWrapper::train(struct svm_problem* problem)
  {
    if (problem != nullptr
       && param_ != nullptr
       && (svm_check_parameter(problem, param_) == nullptr))
    {
      training_set_ = problem;

      if (model_ != nullptr)
      {
#if OPENMS_LIBSVM_VERSION_MAJOR == 2
        svm_destroy_model(model_);
#else
        svm_free_and_destroy_model(&model_);
#endif
        model_ = nullptr;
      }

      if (kernel_type_ == OLIGO)
      {
        if (border_length_ != gauss_table_.size())
        {
          SVMWrapper::calculateGaussTable(border_length_, sigma_, gauss_table_);
        }
        training_problem_ = computeKernelMatrix(problem, problem);
        problem = training_problem_;
      }

      model_ = svm_train(problem, param_);

      return 1;
    }
    else
    {
      if (problem == nullptr)
      {
        cout << "problem is null" << endl;
      }
      if (param_ == nullptr)
      {
        cout << "param_ == null" << endl;
      }
      if (svm_check_parameter(problem, param_) != nullptr)
      {
        cout << "check parameter failed: " << endl
             << svm_check_parameter(problem, param_) << endl;
      }
      cout << "Training error" << endl;
      return 0;
    }
  }

  Int SVMWrapper::train(SVMData& problem)
  {
    if (param_ != nullptr || kernel_type_ != OLIGO)
    {
      training_data_ = problem;

      if (model_ != nullptr)
      {
#if OPENMS_LIBSVM_VERSION_MAJOR == 2
        svm_destroy_model(model_);
#else
        svm_free_and_destroy_model(&model_);
#endif
        model_ = nullptr;
      }

      if (border_length_ != gauss_table_.size())
      {
        SVMWrapper::calculateGaussTable(border_length_, sigma_, gauss_table_);
      }
      training_problem_ = computeKernelMatrix(problem, problem);

      if (svm_check_parameter(training_problem_, param_) == nullptr)
      {
        model_ = svm_train(training_problem_, param_);
        return 1;
      }
    }
    if (training_problem_ == nullptr)
    {
      cout << "problem is null" << endl;
    }
    if (param_ == nullptr)
    {
      cout << "param_ == null" << endl;
    }
    if (svm_check_parameter(training_problem_, param_) != nullptr)
    {
      cout << "check parameter failed" << endl;
    }
    cout << "Training error" << endl;
    return 0;

  }

  void SVMWrapper::saveModel(std::string model_filename) const
  {
    Int  status = 0;

    if (model_ != nullptr)
    {
      status = svm_save_model(model_filename.c_str(), model_);
    }
    else
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, model_filename);
    }
    if (status == -1)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, model_filename);
    }
  }

  void SVMWrapper::loadModel(std::string model_filename)
  {
    TextFile file;
    TextFile::ConstIterator it;
    vector<String> parts;

    if (model_ != nullptr)
    {
#if OPENMS_LIBSVM_VERSION_MAJOR == 2
      svm_destroy_model(model_);
#else
      svm_free_and_destroy_model(&model_);
#endif
      model_ = nullptr;
    }
    model_ = svm_load_model(model_filename.c_str());
    setParameter(SVM_TYPE, svm_get_svm_type(model_));
    file.load(model_filename, true);

    it = StringListUtils::searchPrefix(file.begin(), file.end(), "kernel_type");
    if (it != file.end())
    {
      it->split(' ', parts);
      if (parts[1] == "linear")
      {
        setParameter(KERNEL_TYPE, LINEAR);
      }
      else if (parts[1] == "polynomial")
      {
        setParameter(KERNEL_TYPE, POLY);
      }
      else if (parts[1] == "rbf")
      {
        setParameter(KERNEL_TYPE, RBF);
      }
      else if (parts[1] == "sigmoid")
      {
        setParameter(KERNEL_TYPE, SIGMOID);
      }
      else if (parts[1] == "precomputed")
      {
        setParameter(KERNEL_TYPE, OLIGO);
      }
    }
  }

  void SVMWrapper::predict(struct svm_problem* problem, vector<double>& results)
  {
    results.clear();

    if (model_ == nullptr)
    {
      cout << "Model is null" << endl;
    }
    if (problem == nullptr)
    {
      cout << "problem is null" << endl;
    }
    if (param_->kernel_type == PRECOMPUTED && training_set_ == nullptr)
    {
      cout << "Training set is null and kernel type == PRECOMPUTED" << endl;
    }

    if (model_ != nullptr && problem != nullptr)
    {
      if (kernel_type_ == OLIGO)
      {
        if (training_set_ != nullptr)
        {
          problem = computeKernelMatrix(problem, training_set_);
        }
      }
      results.reserve(problem->l);
      for (Int i = 0; i < problem->l; i++)
      {
        double label = svm_predict(model_, problem->x[i]);
        results.push_back(label);
      }

      if (kernel_type_ == OLIGO)
      {
        LibSVMEncoder::destroyProblem(problem);
      }
    }
  }

  void SVMWrapper::predict(const SVMData& problem, vector<double>& results)
  {
    results.clear();

    if (kernel_type_ == OLIGO)
    {
      if (model_ == nullptr)
      {
        cout << "Model is null" << endl;
      }
      else if (problem.sequences.empty())
      {
        cout << "problem is empty" << endl;
      }
      else if (training_data_.sequences.empty())
      {
        cout << "Training set is empty and kernel type == PRECOMPUTED" << endl;
      }
      else if (model_ != nullptr)
      {
        struct svm_problem* prediction_problem = computeKernelMatrix(problem, training_data_);
        for (Size i = 0; i < problem.sequences.size(); i++)
        {
          double label = svm_predict(model_, prediction_problem->x[i]);
          results.push_back(label);
        }

        LibSVMEncoder::destroyProblem(prediction_problem);
      }
    }
  }

  void SVMWrapper::createRandomPartitions(svm_problem* problem,
                                          Size                                number,
                                          vector<svm_problem*>& problems)
  {
    vector<Size> indices;
    vector<Size>::iterator indices_iterator;

    for (Size i = 0; i < problems.size(); ++i)
    {
      delete problems[i];
    }
    problems.clear();

    if (number == 1)
    {
      problems.push_back(problem);
    }
    else if (number > 1)
    {
      // Creating the particular partition instances
      for (Size i = 0; i < number; i++)
      {
        problems.push_back(new svm_problem());
      }

      // Creating indices
      for (Int  i = 0; i < problem->l; i++)
      {
        indices.push_back(i);
      }
      // Shuffling the indices => random indices
      random_shuffle(indices.begin(), indices.end());

      indices_iterator = indices.begin();

      for (Size partition_index = 0;
           partition_index < number;
           partition_index++)
      {
        Size actual_partition_size = 0;
        // determining the number of elements in this partition
        Size partition_count = (problem->l / number);
        if (problem->l % number > partition_index)
        {
          partition_count++;
        }

        // filling the actual partition with 'partition_count' elements
        while (actual_partition_size < partition_count)
        {
          if (actual_partition_size == 0)
          {
            problems[partition_index]->l = (Int)partition_count;
            problems[partition_index]->x = new svm_node*[partition_count];
            problems[partition_index]->y = new double[partition_count];
          }
          problems[partition_index]->x[actual_partition_size] =
            problem->x[*indices_iterator];
          problems[partition_index]->y[actual_partition_size] =
            problem->y[*indices_iterator];
          ++actual_partition_size;
          ++indices_iterator;
        }
      }
    }
  }

  void SVMWrapper::createRandomPartitions(const SVMData& problem,
                                          Size                                  number,
                                          vector<SVMData>& problems)
  {
    vector<Size> indices;
    vector<Size>::iterator indices_iterator;

    for (Size i = 0; i < problems.size(); ++i)
    {
      problems[i].labels.clear();
      problems[i].sequences.clear();
    }
    problems.clear();

    if (number == 1)
    {
      problems.push_back(problem);
    }
    else if (number > 1)
    {
      // Creating the particular partition instances
      for (Size i = 0; i < number; i++)
      {
        problems.push_back(SVMData());
      }

      // Creating indices
      for (Size  i = 0; i < problem.sequences.size(); i++)
      {
        indices.push_back(i);
      }
      // Shuffling the indices => random indices
      random_shuffle(indices.begin(), indices.end());

      indices_iterator = indices.begin();

      for (Size partition_index = 0;
           partition_index < number;
           partition_index++)
      {
        Size actual_partition_size = 0;
        // determining the number of elements in this partition
        Size partition_count = (problem.sequences.size() / number);
        if (problem.sequences.size() % number > partition_index)
        {
          partition_count++;
        }

        // filling the actual partition with 'partition_count' elements
        while (actual_partition_size < partition_count)
        {
          if (actual_partition_size == 0)
          {
            problems[partition_index].sequences.resize(partition_count, std::vector<std::pair<int, double> >());
            problems[partition_index].labels.resize(partition_count, 0.);
          }
          problems[partition_index].sequences[actual_partition_size] =
            problem.sequences[*indices_iterator];
          problems[partition_index].labels[actual_partition_size] =
            problem.labels[*indices_iterator];
          ++actual_partition_size;
          ++indices_iterator;
        }
      }
    }
  }

  svm_problem* SVMWrapper::mergePartitions(const vector<svm_problem*>& problems, Size except)
  {
    svm_problem* merged_problem = nullptr;

    if (problems.size() == 1 && except == 0)
    {
      return nullptr;
    }

    if (problems.size() > 0)
    {
      int count = 0;

      merged_problem = new svm_problem();
      for (Size i = 0; i < problems.size(); i++)
      {
        if (i != except)
        {
          count += problems[i]->l;
        }
      }
      merged_problem->l = count;
      merged_problem->x = new svm_node*[count];
      merged_problem->y = new double[count];

      Size actual_index = 0;

      for (Size i = 0; i < problems.size(); i++)
      {
        if (i != except)
        {
          for (Int  j = 0; j < problems[i]->l; j++)
          {
            merged_problem->x[actual_index] = problems[i]->x[j];
            merged_problem->y[actual_index] = problems[i]->y[j];
            actual_index++;
          }
        }
      }
    }
    return merged_problem;
  }

  void SVMWrapper::mergePartitions(const vector<SVMData>& problems,
                                   Size                                 except,
                                   SVMData& merged_problem)
  {



    merged_problem.sequences.clear();
    merged_problem.labels.clear();

    if (problems.size() != 1 || except != 0)
    {
      if (problems.size() > 0)
      {
        Size count = 0;
        for (Size i = 0; i < problems.size(); i++)
        {
          if (i != except)
          {
            count += problems[i].labels.size();
          }
        }
        merged_problem.sequences.resize(count, std::vector<std::pair<int, double> >());
        merged_problem.labels.resize(count, 0.);

        Size actual_index = 0;

        for (Size i = 0; i < problems.size(); i++)
        {
          if (i != except)
          {
            for (Size j = 0; j < problems[i].sequences.size(); j++)
            {
              merged_problem.sequences[actual_index] = problems[i].sequences[j];
              merged_problem.labels[actual_index] = problems[i].labels[j];
              actual_index++;
            }
          }
        }
      }
    }
  }

  void SVMWrapper::getLabels(svm_problem* problem, vector<double>& labels)
  {
    labels.clear();

    if (problem != nullptr)
    {
      int count = problem->l;
      for (int i = 0; i < count; i++)
      {
        labels.push_back(problem->y[i]);
      }
    }
  }

  /**
        NEW method:

  */
  double SVMWrapper::performCrossValidation(svm_problem* problem_ul,
                                            const SVMData& problem_l,
                                            const bool                                   is_labeled,
                                            const   map<SVM_parameter_type, double>& start_values_map,
                                            const   map<SVM_parameter_type, double>& step_sizes_map,
                                            const   map<SVM_parameter_type, double>& end_values_map,
                                            Size                                                                             number_of_partitions,
                                            Size                                                                             number_of_runs,
                                            map<SVM_parameter_type, double>& best_parameters,
                                            bool                                                                                 additive_step_sizes,
                                            bool                                                                             output,
                                            String                                                                           performances_file_name,
                                            bool                                                                                 mcc_as_performance_measure)
  {
    map<SVM_parameter_type, double>::const_iterator start_values_iterator;
    vector<pair<double, Size> > combined_parameters;
    combined_parameters.push_back(make_pair(1, 25));

    double cv_quality = 0.0;
    for (Size i = 1; i < gauss_tables_.size(); ++i)
    {
      combined_parameters.push_back(make_pair(1, 25));
    }

    std::vector<double> start_values(start_values_map.size());
    std::vector<double> actual_values(start_values_map.size());
    std::vector<double> step_sizes(start_values_map.size());
    std::vector<double> end_values(start_values_map.size());
    std::vector<SVM_parameter_type> actual_types(start_values_map.size());
    std::vector<double> best_values(start_values_map.size());

    bool found = false; // does a valid grid search cell (with a certain parameter combination) exist?
    Size counter = 0;
    vector<svm_problem*> partitions_ul;
    svm_problem** training_data_ul = nullptr;
    vector<SVMData> partitions_l;
    vector<SVMData> training_data_l;
    double temp_performance = 0;
    vector<double> predicted_labels;
    vector<double> real_labels;
    vector<double> performances;
    Size max_index = 0;
    double max = 0;
    ofstream performances_file;
    ofstream run_performances_file;

    best_parameters.clear();

    if (output)
    {
      performances_file.open(performances_file_name.c_str(), ios_base::out);
      run_performances_file.open((performances_file_name + "_runs.txt").c_str(), ios_base::out);
    }


    // Initializing the grid search variables
    Size actual_index = 0;
    start_values_iterator = start_values_map.begin();
    while (start_values_iterator != start_values_map.end())
    {
      actual_types[actual_index] = start_values_iterator->first;
      actual_values[actual_index] = start_values_iterator->second;
      start_values[actual_index] = start_values_iterator->second;

      map<SVM_parameter_type, double>::const_iterator it = step_sizes_map.find(start_values_iterator->first);
      if (it == step_sizes_map.end())
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No step size given for svm parameter grid search");
      else
        step_sizes[actual_index] = it->second;

      it = end_values_map.find(start_values_iterator->first);
      if (it == end_values_map.end())
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No end value given for svm parameter grid search");
      else
        end_values[actual_index] = it->second;

      ++start_values_iterator;
      ++actual_index;
    }

    // estimate number of training runs:
    Size work_steps(1), work_steps_count(0);
    while (nextGrid_(start_values, step_sizes, end_values, additive_step_sizes, actual_values))
      ++work_steps;
    //reset actual values:
    actual_values = start_values;
    OPENMS_LOG_INFO << "SVM-CrossValidation -- number of grid cells:" << work_steps << "\n";

    work_steps *= number_of_runs * number_of_partitions;
    startProgress(0, work_steps, "SVM-CrossValidation");

    // for every run (each run is identical, except for random partitioning of the data)
    for (Size i = 0; i < number_of_runs; i++)
    {
      for (Size index = 0; index < start_values_map.size(); ++index)
      {
        best_values[index] = 0;
      }
      double max_performance = 0;
      if (is_labeled)
        createRandomPartitions(problem_l, number_of_partitions, partitions_l);
      else
        createRandomPartitions(problem_ul, number_of_partitions, partitions_ul);

      counter = 0;
      found = true;

      if (is_labeled)
        training_data_l.resize(number_of_partitions, SVMData());
      else
        training_data_ul = new svm_problem*[number_of_partitions];
      for (Size j = 0; j < number_of_partitions; j++)
      {
        if (is_labeled)
          SVMWrapper::mergePartitions(partitions_l, j, training_data_l[j]);
        else
          training_data_ul[j] = SVMWrapper::mergePartitions(partitions_ul, j);
      }

      while (found) // do grid search
      {
        // setting svm parameters
        for (Size v = 0; v < start_values_map.size(); ++v)
        {
          // testing whether actual parameters are in the defined range
          if (actual_values[v] > end_values[v])
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RTModel CV parameters are out of range!");
          setParameter(actual_types[v], actual_values[v]);
        }

        temp_performance = 0;

        vector<double>::iterator it_start, it_end;

        // loop over PARTITIONS
        for (Size j = 0; j < number_of_partitions; j++)
        {
          setProgress(work_steps_count++);

          bool success;
          if (is_labeled)
            success = train(training_data_l[j]);
          else
            success = train(training_data_ul[j]);

          if (success)
          {
            if (is_labeled)
            {
              predict(partitions_l[j], predicted_labels);

              it_start = partitions_l[j].labels.begin();
              it_end = partitions_l[j].labels.end();
            }
            else
            {
              predict(partitions_ul[j], predicted_labels);
              getLabels(partitions_ul[j], real_labels);

              it_start = real_labels.begin();
              it_end = real_labels.end();
            }

            if (param_->svm_type == C_SVC || param_->svm_type == NU_SVC)
            {
              if (mcc_as_performance_measure)
              {
                temp_performance +=
                  OpenMS::Math::matthewsCorrelationCoefficient(predicted_labels.begin(), predicted_labels.end(), it_start, it_end);
              }
              else
              {
                temp_performance +=
                  OpenMS::Math::classificationRate(predicted_labels.begin(), predicted_labels.end(), it_start, it_end);
              }
            }
            else if (param_->svm_type == NU_SVR || param_->svm_type == EPSILON_SVR)
            {
              temp_performance +=
                Math::pearsonCorrelationCoefficient(predicted_labels.begin(), predicted_labels.end(), it_start, it_end);
            }

            if (param_->kernel_type == PRECOMPUTED)
            {
              LibSVMEncoder::destroyProblem(training_problem_);
            }

            if (output && j == number_of_partitions - 1)
            {
              performances_file << temp_performance / (j + 1) << " ";
              for (Size k = 0; k < start_values_map.size(); k++)
              {
                switch (actual_types[k])
                {
                case C:
                  performances_file << "C: " << actual_values[k];
                  break;

                case NU:
                  performances_file << "NU: " << actual_values[k];
                  break;

                case DEGREE:
                  performances_file << "DEGREE: " << actual_values[k];
                  break;

                case P:
                  performances_file << "P: " << actual_values[k];
                  break;

                case GAMMA:
                  performances_file << "GAMMA: " << actual_values[k];
                  break;

                case SIGMA:
                  performances_file << "SIGMA: " << actual_values[k];
                  break;

                default:
                  break;
                }
                if (k < (start_values_map.size() - 1))
                {
                  performances_file << " ";
                }
                else
                {
                  performances_file << endl;
                }
              }
            }
          }
          else
          {
            cout << "Training failed" << endl;
          }
        } // ! partitions

        // storing performance for this parameter combination
        temp_performance = temp_performance / number_of_partitions;
        if (temp_performance > max_performance)
        {
          max_performance = temp_performance;
          for (Size index = 0; index < start_values_map.size(); ++index)
          {
            best_values[index] = actual_values[index];
          }
        }

        if (i == 0) // 1st run, append performance for each parameter setting
        {
          performances.push_back(temp_performance);
        }
        else // 2nd+ run, add performance (will be averaged later)
        {
          performances[counter] = performances[counter] + temp_performance;
          counter++;
        }

        // trying to find new parameter combination
        found = nextGrid_(start_values, step_sizes, end_values, additive_step_sizes, actual_values);
      } // ! grid search

      if (!is_labeled)
      {
        for (Size k = 0; k < number_of_partitions; k++)
        {
          delete[] training_data_ul[k]->x;
          delete[] training_data_ul[k]->y;
          delete training_data_ul[k]; // delete individual objects
        }
        delete[] training_data_ul; // delete array of pointers
      }

      // not essential...
      if (output)
      {
        cout << "run finished, time elapsed since start: " << clock()
             << " mean performance is: " << *(max_element(performances.begin(), performances.end())) / (i + 1) << endl
             << "performance of this run is: " << max_performance << " with parameters: ";
        run_performances_file << max_performance << " ";
        for (Size k = 0; k < start_values_map.size(); k++)
        {
          switch (actual_types[k])
          {
          case C:
            cout << "C: " << best_values[k];
            run_performances_file << "C: " << best_values[k];
            break;

          case NU:
            cout << "NU: " << best_values[k];
            run_performances_file << "NU: " << best_values[k];
            break;

          case DEGREE:
            cout << "DEGREE: " << best_values[k];
            run_performances_file << "DEGREE: " << best_values[k];
            break;

          case P:
            cout << "P: " << best_values[k];
            run_performances_file << "P: " << best_values[k];
            break;

          case GAMMA:
            cout << "GAMMA: " << best_values[k];
            run_performances_file << "GAMMA: " << best_values[k];
            break;

          case SIGMA:
            cout << "SIGMA: " << best_values[k];
            run_performances_file << "SIGMA: " << best_values[k];
            break;

          default:
            break;
          }
          if (k < (start_values_map.size() - 1))
          {
            cout << " ";
            run_performances_file << " ";
          }
          else
          {
            cout << endl;
            run_performances_file << endl;
          }
        }
      } // !output

    } // ! number_of_runs "i"

    // Determining the index for the maximum performance
    for (Size i = 0; i < performances.size(); i++)
    {
      if (performances[i] > max)
      {
        max_index = i;
        max = performances[i];
      }
    }

    // Determining the best parameter combination (not necessarily the best single performance, rather average performance of a parameter
    // setting over all runs
    actual_index = 0;
    // resetting actual values to start values (to determine the actual content of the grid cell, we only know its index 'max_index')
    for (start_values_iterator = start_values_map.begin(); start_values_iterator != start_values_map.end(); ++start_values_iterator)
    {
      actual_values[actual_index++] = start_values_iterator->second;
    }

    if (max_index == 0)
    {
      // do nothing, the initial parameters were the best
    }
    else
    {
      performances_file << "Best parameter combination *********************" << endl;
      counter = 1;
      found = true;
      while (found)
      {

        if (counter <= max_index)
          found = nextGrid_(start_values, step_sizes, end_values, additive_step_sizes, actual_values);

        performances_file << performances[counter]  / number_of_runs << ": ";
        for (Size k = 0; k < start_values_map.size(); k++)
        {
          switch (actual_types[k])
          {
          case C:
            performances_file << "C: " << actual_values[k];
            break;

          case NU:
            performances_file << "NU: " << actual_values[k];
            break;

          case DEGREE:
            performances_file << "DEGREE: " << actual_values[k];
            break;

          case P:
            performances_file << "P: " << actual_values[k];
            break;

          case GAMMA:
            performances_file << "GAMMA: " << actual_values[k];
            break;

          case SIGMA:
            performances_file << "SIGMA: " << actual_values[k];
            break;

          default:
            break;
          }
          if (k < (start_values_map.size() - 1))
          {
            performances_file << ", ";
          }
          else
          {
            performances_file << endl;
          }
        }
        // best parameter combination found
        if (counter == max_index)
        {
          found = false;
        }
        counter++;
      }
      performances_file << "Best parameter combination ended****************" << endl;
    }

    if (output)
    {
      performances_file.flush();
      run_performances_file.flush();
      performances_file.close();
      run_performances_file.close();
    }

    for (actual_index = 0; actual_index < start_values_map.size(); actual_index++)
    {
      best_parameters.insert(make_pair(actual_types[actual_index], actual_values[actual_index]));
    }
    cv_quality = performances[max_index] / number_of_runs;

    actual_index = 0;
    while (actual_index < start_values_map.size())
    {
      setParameter(actual_types[actual_index], actual_values[actual_index]);
      actual_index++;
    }

    return cv_quality;

  }

  bool SVMWrapper::nextGrid_(const std::vector<double>& start_values,
                             const std::vector<double>& step_sizes,
                             const std::vector<double>& end_values,
                             const bool additive_step_sizes,
                             std::vector<double>& actual_values)
  {
    bool found = false;
    Size actual_index = 0;
    double precision = 0.0001;

    while (actual_index < start_values.size() && !found)
    {
      double new_value;
      if (additive_step_sizes)
        new_value = actual_values[actual_index] + step_sizes[actual_index];
      else
        new_value = actual_values[actual_index] * step_sizes[actual_index];

      if (new_value   <= end_values[actual_index] + precision)
      {
        found = true;
        actual_values[actual_index] = new_value;
      }
      else
      {
        actual_values[actual_index] = start_values[actual_index];
      }

      actual_index++;
    }

    return found;
  }

  void SVMWrapper::predict(const std::vector<svm_node*>& vectors, vector<double>& results)
  {
    results.clear();

    if (model_ != nullptr)
    {
      for (Size i = 0; i < vectors.size(); i++)
      {
        double label = svm_predict(model_, vectors[i]);
        results.push_back(label);
      }
    }
  }

  double SVMWrapper::getSVRProbability()
  {
    if (model_ != nullptr)
    {
      return svm_get_svr_probability(model_);
    }
    else
    {
      return 0;
    }
  }

  void SVMWrapper::getSVCProbabilities(struct svm_problem* problem,
                                       vector<double>& probabilities,
                                       vector<double>& prediction_labels)
  {
    vector<double> temp_prob_estimates;
    temp_prob_estimates.resize(2, -1);
    vector<int> labels;
    labels.push_back(-1);
    labels.push_back(1);

    svm_get_labels(model_, &(labels[0]));

    probabilities.clear();
    prediction_labels.clear();

    if (model_ != nullptr)
    {
      if (kernel_type_ == OLIGO)
      {
        if (training_set_ != nullptr)
        {
          problem = computeKernelMatrix(problem, training_set_);
        }
      }
      for (int i = 0; i < problem->l; ++i)
      {
        double temp_label = svm_predict_probability(model_, problem->x[i], &(temp_prob_estimates[0]));
        prediction_labels.push_back(temp_label);
        if (labels[0] >= 0)
        {
          probabilities.push_back(temp_prob_estimates[0]);
        }
        else
        {
          probabilities.push_back(1 - temp_prob_estimates[0]);
        }
      }
      if (kernel_type_ == OLIGO)
      {
        LibSVMEncoder::destroyProblem(problem);
      }
    }
  }

  void SVMWrapper::initParameters_()
  {
    model_ = nullptr;

    param_->svm_type = NU_SVR;
    param_->kernel_type = PRECOMPUTED;
    param_->degree = 1; // for poly
    param_->gamma = 1.0; // for poly/rbf/sigmoid
    param_->coef0 = 0; // for poly/sigmoid
    param_->cache_size = 300; // in MB
    param_->eps = 0.001; // stopping criterion
    param_->C = 1; // for C_SVC, EPSILON_SVR, and NU_SVR
    param_->nu = 0.5; // for NU_SVC, ONE_CLASS, and NU_SVR
    param_->p = 0.1; // for EPSILON_SVR
    param_->shrinking = 0; // use the shrinking heuristics
    param_->probability = 0; // do not build probability model during training

    param_->nr_weight = 0; // for C_SVC
    param_->weight_label = nullptr; // for C_SVC
    param_->weight = nullptr; // for C_SVC

    // silence libsvm
    svm_set_print_string_function(&printToVoid_);
  }

  void SVMWrapper::printToVoid_(const char* /*s*/)
  {
    return;
  }

  void SVMWrapper::setWeights(const vector<Int>& weight_labels, const vector<double>& weights)
  {
    if (weight_labels.size() == weights.size() && weights.size() > 0)
    {
      param_->nr_weight = (Int)weights.size();
      param_->weight_label = new Int[weights.size()];
      param_->weight = new double[weights.size()];
      for (Size i = 0; i < weights.size(); ++i)
      {
        param_->weight_label[i] = weight_labels[i];
        param_->weight[i] = weights[i];
      }
    }
  }

  double SVMWrapper::kernelOligo(const vector<pair<int, double> >& x,
                                 const vector<pair<int, double> >& y,
                                 const vector<double>& gauss_table,
                                 int                                                                 max_distance)
  {
    double kernel = 0;
    Size i1     = 0;
    Size i2     = 0;
    Size c1     = 0;
    Size x_size = x.size();
    Size y_size = y.size();

    while (i1 < x_size && i2 < y_size)
    {
      if (x[i1].second == y[i2].second)
      {
        if (max_distance < 0
           || (abs(x[i1].first - y[i2].first)) <= max_distance)
        {
          kernel += gauss_table.at(abs((x[i1].first - y[i2].first)));
          if (i1 < x_size - 1 && x[i1].second == x[i1 + 1].second)
          {
            i1++;
            c1++;
          }
          else if (i2 < y_size - 1 && y[i2].second == y[i2 + 1].second)
          {
            i2++;
            i1 -= c1;
            c1 = 0;
          }
          else
          {
            i1++;
            i2++;
          }
        }
        else
        {
          if (x[i1].first < y[i2].first)
          {
            if (i1 < x_size - 1 && x[i1].second == x[i1 + 1].second)
            {
              i1++;
            }
            else if (i2 < y_size - 1 && y[i2].second == y[i2 + 1].second)
            {
              while (i2 < y_size - 1 && y[i2].second == y[i2 + 1].second)
              {
                ++i2;
              }
              ++i1;
              c1 = 0;
            }
            else
            {
              i1++;
              i2++;
              c1 = 0;
            }
          }
          else
          {
            i2++;
            i1 -= c1;
            c1 = 0;
          }
        }
      }
      else
      {
        if (x[i1].second < y[i2].second)
        {
          i1++;
        }
        else
        {
          i2++;
        }
        c1 = 0;
      }
    }
    return kernel;
  }

  double SVMWrapper::kernelOligo(const svm_node* x,
                                 const svm_node* y,
                                 const vector<double>& gauss_table,
                                 double                                 sigma_square,
                                 Size                                           max_distance)
  {
    double kernel = 0;
    Int    i1     = 0;
    Int    i2     = 0;
    Int    c1     = 0;

    while (x[i1].index != -1
          && y[i2].index != -1)
    {
      if (x[i1].index == y[i2].index)
      {
        if (((Size) abs(x[i1].value - y[i2].value)) <= max_distance)
        {
          if (sigma_square == 0)
          {
            try
            {
              kernel += gauss_table.at(abs((int)(x[i1].value
                                                 - y[i2].value)));
            }
            catch (...)
            {
              cout << "Tried to access " << x[i1].value << " - " << y[i2].value << endl;
            }
          }
          else
          {
            kernel += exp(-1 * (x[i1].value - y[i2].value) * (x[i1].value - y[i2].value) / (4 * sigma_square));

          }
          if (x[i1].index == x[i1 + 1].index)
          {
            i1++;
            c1++;
          }
          else if (y[i2].index == y[i2 + 1].index)
          {
            i2++;
            i1 -= c1;
            c1 = 0;
          }
          else
          {
            i1++;
            i2++;
          }
        }
        else
        {
          if (x[i1].value < y[i2].value)
          {
            if (x[i1].index == x[i1 + 1].index)
            {
              i1++;
            }
            else if (y[i2].index == y[i2 + 1].index)
            {
              i2++;
              i1 -= c1;
              c1 = 0;
            }
            else
            {
              i1++;
              i2++;
            }
          }
          else
          {
            i2++;
            i1 -= c1;
            c1 = 0;
          }
        }
      }
      else
      {
        if (x[i1].index < y[i2].index)
        {
          i1++;
        }
        else
        {
          i2++;
        }
        c1 = 0;
      }
    }
    return kernel;
  }

  svm_problem* SVMWrapper::computeKernelMatrix(svm_problem* problem1, svm_problem* problem2)
  {
    double temp = 0;
    svm_problem* kernel_matrix;

    if (problem1 == nullptr || problem2 == nullptr)
    {
      return nullptr;
    }
    UInt number_of_sequences = problem1->l;
    kernel_matrix = new svm_problem;
    kernel_matrix->l = number_of_sequences;
    kernel_matrix->x = new svm_node*[number_of_sequences];
    kernel_matrix->y = new double[number_of_sequences];

    for (Size i = 0; i < number_of_sequences; i++)
    {
      kernel_matrix->x[i] = new svm_node[problem2->l + 2];
      kernel_matrix->x[i][0].index = 0;
      kernel_matrix->x[i][0].value = i + 1;
      kernel_matrix->y[i] = problem1->y[i];
      kernel_matrix->x[i][problem2->l + 1].index = -1;
    }

    if (problem1 == problem2)
    {
      for (Size i = 0; i < number_of_sequences; i++)
      {
        for (Size j = i; j < number_of_sequences; j++)
        {
          temp = SVMWrapper::kernelOligo(problem1->x[i], problem2->x[j], gauss_table_);
          kernel_matrix->x[i][j + 1].index = (Int)j + 1;
          kernel_matrix->x[i][j + 1].value = temp;
          kernel_matrix->x[j][i + 1].index = (Int)i + 1;
          kernel_matrix->x[j][i + 1].value = temp;
        }
      }
    }
    else
    {
      for (Size i = 0; i < number_of_sequences; i++)
      {
        for (Size j = 0; j < (Size) problem2->l; j++)
        {
          temp = SVMWrapper::kernelOligo(problem1->x[i], problem2->x[j], gauss_table_);

          kernel_matrix->x[i][j + 1].index = (Int)j + 1;
          kernel_matrix->x[i][j + 1].value = temp;
        }
      }
    }
    return kernel_matrix;
  }

  svm_problem* SVMWrapper::computeKernelMatrix(const SVMData& problem1, const SVMData& problem2)
  {
    double temp = 0;
    svm_problem* kernel_matrix;

    if (problem1.labels.empty() || problem2.labels.empty())
    {
      return nullptr;
    }

    if (problem1.labels.size() != problem1.sequences.size()
       || problem2.labels.size() != problem2.sequences.size())
    {
      return nullptr;
    }

    Size number_of_sequences = problem1.labels.size();
    kernel_matrix = new svm_problem;
    kernel_matrix->l = (int) number_of_sequences;
    kernel_matrix->x = new svm_node*[number_of_sequences];
    kernel_matrix->y = new double[number_of_sequences];

    for (Size i = 0; i < number_of_sequences; i++)
    {
      kernel_matrix->x[i] = new svm_node[problem2.labels.size() + 2];
      kernel_matrix->x[i][0].index = 0;
      kernel_matrix->x[i][0].value = i + 1;
      kernel_matrix->y[i] = problem1.labels[i];
      kernel_matrix->x[i][problem2.labels.size() + 1].index = -1;
    }

    if (&problem1 == &problem2)
    {
      for (Size i = 0; i < number_of_sequences; i++)
      {
        for (Size j = i; j < number_of_sequences; j++)
        {
          temp = SVMWrapper::kernelOligo(problem1.sequences[i], problem2.sequences[j], gauss_table_);
          kernel_matrix->x[i][j + 1].index = int(j) + 1;
          kernel_matrix->x[i][j + 1].value = temp;
          kernel_matrix->x[j][i + 1].index = int(i) + 1;
          kernel_matrix->x[j][i + 1].value = temp;
        }
      }
    }
    else
    {
      for (Size i = 0; i < number_of_sequences; i++)
      {
        for (Size j = 0; j < problem2.labels.size(); j++)
        {
          temp = SVMWrapper::kernelOligo(problem1.sequences[i], problem2.sequences[j], gauss_table_);

          kernel_matrix->x[i][j + 1].index = int(j) + 1;
          kernel_matrix->x[i][j + 1].value = temp;
        }
      }
    }
    return kernel_matrix;
  }

  void SVMWrapper::getSignificanceBorders(svm_problem* data,
                                          pair<double, double>& sigmas,
                                          double confidence,
                                          Size number_of_runs,
                                          Size number_of_partitions,
                                          double step_size,
                                          Size max_iterations)
  {
    vector<pair<double, double> > points;
    vector<double>                                      differences;
    vector<svm_problem*>                                partitions;
    svm_problem* training_data;
    vector<double>                                      predicted_labels;
    vector<double>                                      real_labels;
    Size                                                                    counter = 0;
    Size                                                                    target = 0;
    ofstream                                                            file("points.txt");
    double                                                      mean;
    double                                                      intercept = 0;
    double                                                      slope = 0;
    double                                                      step_size1 = 0.;
    double                                                      step_size2 = step_size;
    double                                                      maximum = 0.;
    double                                                      minimum = 0.;


    // creation of points (measured rt, predicted rt)
    for (Size i = 0; i < number_of_runs; ++i)
    {
      createRandomPartitions(data, number_of_partitions, partitions);

      for (Size j = 0; j < number_of_partitions; ++j)
      {
        training_data = SVMWrapper::mergePartitions(partitions, j);
        if (train(training_data))
        {
          predict(partitions[j], predicted_labels);
          getLabels(partitions[j], real_labels);
          vector<double>::iterator pred_it = predicted_labels.begin();
          vector<double>::iterator real_it = real_labels.begin();
          while (pred_it != predicted_labels.end()
                && real_it != real_labels.end())
          {
            points.push_back(make_pair(*real_it, *pred_it));
            differences.push_back(abs(*real_it - *pred_it));
            file << *real_it << " " << *pred_it << endl;
            ++pred_it;
            ++real_it;
          }
        }
      }
    }
    file << flush;

    // trying to find the two line parameters
    target = (Size) Math::round(confidence * points.size());

    mean = accumulate(differences.begin(), differences.end(), 0.0) / differences.size();
    intercept = mean;
    slope = 1;

    step_size1 = (maximum - minimum) * step_size;
    while (target > getNumberOfEnclosedPoints_(intercept, slope, points) && counter < max_iterations)
    {

      cout << "intercept: " << intercept << ", slope: " << slope << " shape contains "
           << ((getNumberOfEnclosedPoints_(intercept, slope, points) / ((double)points.size())) * 100)
           << " % of points" << endl;

      intercept += step_size1;
      slope += step_size2;
      ++counter;
    }
    sigmas.first = intercept;
    sigmas.second = slope;

    cout << "intercept: " << intercept << ", slope: " << slope << " shape contains "
         << ((getNumberOfEnclosedPoints_(intercept, slope, points) / ((double)points.size())) * 100)
         << " % of points" << endl;
  }

  void SVMWrapper::getSignificanceBorders(const SVMData& data,
                                          pair<double, double>& sigmas,
                                          double confidence,
                                          Size number_of_runs,
                                          Size number_of_partitions,
                                          double step_size,
                                          Size max_iterations)
  {
    vector<pair<double, double> > points;
    vector<double>                                      differences;
    vector<SVMData>                                         partitions;
    SVMData                                                             training_data;
    vector<double>                                      predicted_labels;
    Size                                                                    counter = 0;
    Size                                                                    target = 0;
    ofstream                                                            file("points.txt");
    double                                                      mean;
    double                                                      intercept = 0;
    double                                                      slope = 0;
    double                                                      step_size1 = 0.;
    double                                                      step_size2 = step_size;
    double                                                      maximum = 0.;
    double                                                      minimum = 0.;



    // creation of points (measured rt, predicted rt)
    for (Size i = 0; i < number_of_runs; ++i)
    {
      createRandomPartitions(data, number_of_partitions, partitions);

      for (Size j = 0; j < number_of_partitions; ++j)
      {
        SVMWrapper::mergePartitions(partitions, j, training_data);
        if (train(training_data))
        {
          predict(partitions[j], predicted_labels);
          vector<double>::iterator pred_it = predicted_labels.begin();
          vector<double>::iterator real_it = partitions[j].labels.begin();
          while (pred_it != predicted_labels.end()
                && real_it != partitions[j].labels.end())
          {
            points.push_back(make_pair(*real_it, *pred_it));
            differences.push_back(abs(*real_it - *pred_it));
            file << *real_it << " " << *pred_it << endl;

            if (*real_it < minimum)
            {
              minimum = *real_it;
            }
            if (*real_it > maximum)
            {
              maximum = *real_it;
            }
            ++pred_it;
            ++real_it;
          }
        }
      }
    }
    file << flush;

    // trying to find the two line parameters
    target = (Size) Math::round(confidence * points.size());

    mean = accumulate(differences.begin(), differences.end(), 0.0) / differences.size();
    intercept = mean;
    slope = 1;

    step_size1 = (maximum - minimum) * step_size;
    while (target > getNumberOfEnclosedPoints_(intercept, slope, points) && counter < max_iterations)
    {

      cout << "intercept: " << intercept << ", slope: " << slope << " shape contains "
           << ((getNumberOfEnclosedPoints_(intercept, slope, points) / ((double)points.size())) * 100)
           << " % of points" << endl;

      intercept += step_size1;
      slope += step_size2;
      ++counter;
    }
    sigmas.first = intercept;
    sigmas.second = slope;

    cout << "intercept: " << intercept << ", slope: " << slope << " shape contains "
         << ((getNumberOfEnclosedPoints_(intercept, slope, points) / ((double)points.size())) * 100)
         << " % of points" << endl;
  }

  Size SVMWrapper::getNumberOfEnclosedPoints_(double intercept,
                                              double slope,
                                              const vector<pair<double, double> >& points)
  {
    Size counter = 0;

    for (vector<pair<double, double> >::const_iterator it = points.begin();
         it != points.end();
         ++it)
    {
      double upper_bound = intercept + it->first * slope;
      double lower_bound = -intercept + it->first * (1 / slope);

      if (it->second >= lower_bound
         && it->first <= upper_bound)
      {
        ++counter;
      }
    }
    return counter;
  }

  double SVMWrapper::getPValue(double intercept,
                               double slope,
                               pair<double, double> point)
  {
    double center = point.first;
    double distance = std::abs(point.second - center);

    double actual_sigma = ((intercept + point.first * slope) - point.first) / 2;

    actual_sigma = (actual_sigma < 0) ? (-1 * actual_sigma) : actual_sigma;
    double sd_units = distance / actual_sigma;

    // getting only the inner part of the area [-1,1]
    boost::math::normal dist(0, 1);

    return (boost::math::cdf(dist, sd_units) - 0.5) * 2;
  }

  void SVMWrapper::getDecisionValues(svm_problem* data, vector<double>& decision_values)
  {
    double temp_value;

    decision_values.clear();
    if (model_ != nullptr)
    {
      if (param_->svm_type == NU_SVR || param_->svm_type == EPSILON_SVR)
      {
        predict(data, decision_values);
      }
      else if (svm_get_nr_class(model_) == 2)
      {
        bool first_label_positive = false;
        std::vector<int> labels;
        labels.resize(svm_get_nr_class(model_));
        svm_get_labels(model_, &(labels[0]));
        if (labels[0] == 1)
        {
          first_label_positive = true;
        }

        if (kernel_type_ == OLIGO)
        {
          if (training_set_ != nullptr)
          {
            data = computeKernelMatrix(data, training_set_);
          }
        }
        for (Int  i = 0; i < data->l; ++i)
        {
          temp_value = 0;
          svm_predict_values(model_, data->x[i], &temp_value);
          if (first_label_positive)
          {
            decision_values.push_back(temp_value);
          }
          else
          {
            decision_values.push_back(-temp_value);
          }
        }
        if (kernel_type_ == OLIGO)
        {
          LibSVMEncoder::destroyProblem(data);
        }
      }
    }
  }

  void SVMWrapper::scaleData(svm_problem* data, Int  max_scale_value)
  {
    vector<double> max_values;
    vector<double> min_values;
    vector<double> sums;
    Int  max_index = 0;
    Int  j = 0;

    for (Int  i = 0; i < data->l; ++i)
    {
      j = 0;
      while (data->x[i][j].index != -1)
      {
        if (data->x[i][j].index > max_index)
        {
          max_index = data->x[i][j].index;
        }
        ++j;
      }
    }

    max_values.resize(max_index, 0);
    min_values.resize(max_index, 0);
    sums.resize(max_index, 0);

    for (Int  i = 0; i < data->l; ++i)
    {
      j = 0;
      while (data->x[i][j].index != -1)
      {
        if (data->x[i][j].value > max_values.at(data->x[i][j].index - 1))
        {
          max_values.at(data->x[i][j].index - 1) = data->x[i][j].value;
        }
        sums.at(data->x[i][j].index - 1) = sums.at(data->x[i][j].index - 1) + data->x[i][j].value;
        if (data->x[i][j].value < min_values.at(data->x[i][j].index - 1))
        {
          min_values.at(data->x[i][j].index - 1) = data->x[i][j].value;
        }

        ++j;
      }
    }
    for (Int  i = 0; i < data->l; ++i)
    {
      j = 0;
      while (data->x[i][j].index != -1)
      {
        if (max_scale_value == -1)
        {
          data->x[i][j].value = 2 * (data->x[i][j].value - min_values.at(data->x[i][j].index - 1))
                                / (max_values.at(data->x[i][j].index - 1) - min_values.at(data->x[i][j].index - 1)) - 1;
        }
        else
        {
          data->x[i][j].value = max_scale_value * (data->x[i][j].value - min_values.at(data->x[i][j].index - 1))
                                / (max_values.at(data->x[i][j].index - 1) - min_values.at(data->x[i][j].index - 1));
        }
        ++j;
      }
    }
  }

  void SVMWrapper::calculateGaussTable(Size border_length,
                                       double sigma,
                                       vector<double>& gauss_table)
  {
    if (border_length != gauss_table.size())
    {
      gauss_table.resize(border_length, 0);
    }
    gauss_table[0] = 1;
    for (Size i = 1; i < border_length; ++i)
    {
      gauss_table[i] = exp((-1 / 4.0 /
                            (sigma * sigma)) *
                           (i * i));
    }
  }

} // namespace OpenMS
