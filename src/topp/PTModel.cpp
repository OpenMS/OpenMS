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
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_PTModel PTModel

    @brief Used to train a model for the prediction of proteotypic peptides.

    The input consists of two files: One file contains the positive examples (the peptides which
    are proteotypic) and the other contains the negative examples (the nonproteotypic peptides).

    Parts of this model has been described in the publication

    Ole Schulz-Trieglaff, Nico Pfeifer, Clemens Gr&ouml;pl, Oliver Kohlbacher and Knut Reinert
    LC-MSsim - a simulation software for Liquid Chromatography Mass Spectrometry data
    BMC Bioinformatics 2008, 9:423.

    There are a number of parameters which can be changed for the svm (specified in the ini file):
    <ul>
        <li>
            kernel_type: the kernel function (e.g., POLY for the
                polynomial kernel, LINEAR for the linear kernel or RBF for the gaussian kernel); we recommend
                SVMWrapper::OLIGO for our paired oligo-border kernel (POBK)
        </li>
        <li>
            border_length: border length for the POBK
        </li>
        <li>
            k_mer_length: length of the signals considered in the POBK
        </li>
        <li>
            sigma: the amount of positional smoothing for the POBK
        </li>
        <li>
            degree: the degree parameter for the polynomial kernel
        </li>
        <li>
            c: the penalty parameter of the svm
        </li>
        <li>
            nu: the nu parameter for nu-SVC
        </li>
    </ul>

    The last five parameters (sigma, degree, c, nu and p)
    are used in a cross validation (CV) to find the best parameters according to the
    training set. Thus, you have to specify the start value of a
    parameter, the step size in which the parameters should be increased
    and a final value for the particular parameter such that the tested
    parameter is never bigger than the given final value. If you want
    to perform a cross validation, for example, for the parameter c, you
    have to specify <b>c_start</b>, <b>c_step_size</b> and <b>c_stop</b>
    in the ini file. Let's say you want to perform a CV for c from 0.1 to 2
    with step size 0.1. Open up your ini-file with INIFileEditor and modify the fields
    c_start, c_step_size, and c_stop accordingly.

    If the CV should test additional parameters in a certain range
    you just include them analogously to the example above.
    Furthermore, you can specify the number of partitions for the CV with
    <b>number_of_partitions</b> in the ini file and the number of runs
    with <b>number_of_runs</b>.

    <br>
    Consequently you have two choices to use this application:

    <ol>
        <li>
            Set the parameters of the svm: The PTModel application will train
            the svm with the training data and store the svm model.
        </li>
        <li>
            Give a range of parameters for which a CV should be performed:
            The PTModel application will perform a CV to find the best
            parameter combination in the given range and afterwards train
            the svm with the best parameters and the whole training data.
            Then the model is stored.
        </li>
    </ol>

    <br>
    The model can be used in @ref TOPP_PTPredict, to predict the likelihood
    for peptides to be proteotypic.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_PTModel.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_PTModel.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPPTModel :
  public TOPPBase
{
public:
  TOPPPTModel() :
    TOPPBase("PTModel", "Trains a model for the prediction of proteotypic peptides from a training set.")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_positive", "<file>", "", "input file with positive examples");
    setValidFormats_("in_positive", ListUtils::create<String>("idXML"));
    registerInputFile_("in_negative", "<file>", "", "input file with negative examples");
    setValidFormats_("in_negative", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file: the model in libsvm format");
    setValidFormats_("out", ListUtils::create<String>("txt"));
    registerOutputFile_("out_oligo_params", "<file>", "", "output file with additional model parameters when using the OLIGO kernel", false);
    setValidFormats_("out_oligo_params", ListUtils::create<String>("paramXML"));
    registerOutputFile_("out_oligo_trainset", "<file>", "", "output file with the used training dataset when using the OLIGO kernel", false);
    setValidFormats_("out_oligo_trainset", ListUtils::create<String>("txt"));
    registerDoubleOption_("c", "<float>", 1, "the penalty parameter of the svm", false);
    registerStringOption_("svm_type", "<type>", "C_SVC", "the type of the svm (NU_SVC or C_SVC)", false);
    setValidStrings_("svm_type", ListUtils::create<String>("NU_SVC,C_SVC"));
    registerDoubleOption_("nu", "<float>", 0.5, "the nu parameter [0..1] of the svm (for nu-SVR)", false);
    setMinFloat_("nu", 0);
    setMaxFloat_("nu", 1);
    registerStringOption_("kernel_type", "<type>", "OLIGO", "the kernel type of the svm", false);
    setValidStrings_("kernel_type", ListUtils::create<String>("LINEAR,RBF,POLY,OLIGO"));
    registerIntOption_("degree", "<int>", 1, "the degree parameter of the kernel function of the svm (POLY kernel)", false);
    setMinInt_("degree", 1);
    registerIntOption_("border_length", "<int>", 22, "length of the POBK", false);
    setMinInt_("border_length", 1);
    registerIntOption_("k_mer_length", "<int>", 1, "k_mer length of the POBK", false);
    setMinInt_("k_mer_length", 1);
    registerDoubleOption_("sigma", "<float>", 5, "sigma of the POBK", false);
    registerIntOption_("max_positive_count", "<int>", 1000, "quantity of positive samples for training (randomly chosen if smaller than available quantity)", false);
    setMinInt_("max_positive_count", 1);
    registerIntOption_("max_negative_count", "<int>", 1000, "quantity of positive samples for training (randomly chosen if smaller than available quantity)", false);
    setMinInt_("max_negative_count", 1);
    registerFlag_("redundant", "if the input sets are redundant and the redundant peptides should occur more than once in the training set, this flag has to be set");
    registerFlag_("additive_cv", "if the step sizes should be interpreted additively (otherwise the actual value is multiplied with the step size to get the new value");

    addEmptyLine_();
    registerTOPPSubsection_("cv", "Parameters for the grid search / cross validation:");
    registerFlag_("cv:skip_cv", "Has to be set if the cv should be skipped and the model should just be trained with the specified parameters.");
    registerIntOption_("cv:number_of_runs", "<int>", 10, "number of runs for the CV", false);
    setMinInt_("cv:number_of_runs", 1);
    registerIntOption_("cv:number_of_partitions", "<int>", 10, "number of CV partitions", false);
    setMinInt_("cv:number_of_partitions", 2);
    registerIntOption_("cv:degree_start", "<int>", 1, "starting point of degree", false);
    setMinInt_("cv:degree_start", 1);
    registerIntOption_("cv:degree_step_size", "<int>", 2, "step size point of degree", false);
    registerIntOption_("cv:degree_stop", "<int>", 4, "stopping point of degree", false);
    registerDoubleOption_("cv:c_start", "<float>", 1, "starting point of c", false);
    registerDoubleOption_("cv:c_step_size", "<float>", 100, "step size of c", false);
    registerDoubleOption_("cv:c_stop", "<float>", 1000, "stopping point of c", false);
    registerDoubleOption_("cv:nu_start", "<float>", 0.1, "starting point of nu", false);
    setMinFloat_("cv:nu_start", 0);
    setMaxFloat_("cv:nu_start", 1);
    registerDoubleOption_("cv:nu_step_size", "<float>", 1.3, "step size of nu", false);
    registerDoubleOption_("cv:nu_stop", "<float>", 0.9, "stopping point of nu", false);
    setMinFloat_("cv:nu_stop", 0);
    setMaxFloat_("cv:nu_stop", 1);
    registerDoubleOption_("cv:sigma_start", "<float>", 1, "starting point of sigma", false);
    registerDoubleOption_("cv:sigma_step_size", "<float>", 1.3, "step size of sigma", false);
    registerDoubleOption_("cv:sigma_stop", "<float>", 15, "stopping point of sigma", false);
  }

  ExitCodes main_(Int, const char**) override
  {
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    vector<ProteinIdentification> protein_identifications_negative;
    vector<PeptideIdentification> identifications_negative;
    vector<String> training_peptides;
    vector<double> training_labels;
    PeptideHit temp_peptide_hit;
    SVMWrapper svm;
    LibSVMEncoder encoder;
    svm_problem* encoded_training_sample = nullptr;
    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
    map<SVMWrapper::SVM_parameter_type, double> start_values;
    map<SVMWrapper::SVM_parameter_type, double> step_sizes;
    map<SVMWrapper::SVM_parameter_type, double> end_values;
    double sigma_start = 0;
    double sigma_step_size = 0;
    double sigma_stop = 0;
    UInt number_of_partitions = 0;
    UInt number_of_runs = 0;
    map<SVMWrapper::SVM_parameter_type, double> optimized_parameters;
    map<SVMWrapper::SVM_parameter_type, double>::iterator parameters_iterator;
    bool additive_cv = true;
    Param additional_parameters;
    Int temp_type = POLY;
    String debug_string = "";
    double sigma = 0.1;
    UInt k_mer_length = 1;
    Int border_length = 0;
    bool non_redundant = false;
    bool skip_cv = getFlag_("cv:skip_cv");

    svm.setParameter(SVMWrapper::PROBABILITY, 1);
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_positives = getStringOption_("in_positive");
    String inputfile_negatives = getStringOption_("in_negative");
    String temp_string = "";

    String outputfile_name = getStringOption_("out");

    UInt max_positive_count = getIntOption_("max_positive_count");
    UInt max_negative_count = getIntOption_("max_negative_count");

    //SVM type
    String type = getStringOption_("svm_type");
    if (type == "NU_SVC")
    {
      svm.setParameter(SVMWrapper::SVM_TYPE, NU_SVC);
    }
    else if (type == "C_SVC")
    {
      svm.setParameter(SVMWrapper::SVM_TYPE, C_SVC);
    }
    else
    {
      writeLog_("Illegal svm type given. Svm type has to be either "
                + String("NU_SVC or C_SVC. Aborting!"));
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    //Kernel type
    type = getStringOption_("kernel_type");
    if (type == "POLY")
    {
      svm.setParameter(SVMWrapper::KERNEL_TYPE, POLY);
      temp_type = POLY;
    }
    else if (type == "LINEAR")
    {
      svm.setParameter(SVMWrapper::KERNEL_TYPE, LINEAR);
      temp_type = LINEAR;
    }
    else if (type == "RBF")
    {
      svm.setParameter(SVMWrapper::KERNEL_TYPE, RBF);
      temp_type = RBF;
    }
    else if (type == "OLIGO")
    {
      svm.setParameter(SVMWrapper::KERNEL_TYPE, SVMWrapper::OLIGO);
      temp_type = SVMWrapper::OLIGO;
    }
    else if (type == "SIGMOID")
    {
      svm.setParameter(SVMWrapper::KERNEL_TYPE, SIGMOID);
      temp_type = SIGMOID;
    }
    else
    {
      writeLog_("Unknown kernel type given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    //parameters
    svm.setParameter(SVMWrapper::C, getDoubleOption_("c"));
    svm.setParameter(SVMWrapper::DEGREE, getIntOption_("degree"));
    if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVC)
    {
      svm.setParameter(SVMWrapper::NU, getDoubleOption_("nu"));
    }

    //grid search parameters
    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == POLY)
    {
      svm.setParameter(SVMWrapper::DEGREE, getIntOption_("degree"));
      if (!skip_cv)
      {
        double degree_start = getIntOption_("cv:degree_start");
        double degree_step_size = getIntOption_("cv:degree_step_size");
        if (!additive_cv && degree_step_size <= 1)
        {
          writeLog_("Step size of degree <= 1 and additive_cv is false. Aborting!");
          return ILLEGAL_PARAMETERS;
        }
        double degree_stop = getIntOption_("cv:degree_stop");

        start_values.insert(make_pair(SVMWrapper::DEGREE, degree_start));
        step_sizes.insert(make_pair(SVMWrapper::DEGREE, degree_step_size));
        end_values.insert(make_pair(SVMWrapper::DEGREE, degree_stop));
      }
    }

    if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == C_SVC && !skip_cv)
    {
      double c_start = getDoubleOption_("cv:c_start");
      double c_step_size = getDoubleOption_("cv:c_step_size");
      if (!additive_cv && c_step_size <= 1)
      {
        writeLog_("Step size of c <= 1 and additive_cv is false. Aborting!");
        return ILLEGAL_PARAMETERS;
      }
      double c_stop = getDoubleOption_("cv:c_stop");

      start_values.insert(make_pair(SVMWrapper::C, c_start));
      step_sizes.insert(make_pair(SVMWrapper::C, c_step_size));
      end_values.insert(make_pair(SVMWrapper::C, c_stop));
    }

    if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVC && !skip_cv)
    {
      double nu_start = getDoubleOption_("cv:nu_start");
      double nu_step_size = getDoubleOption_("cv:nu_step_size");
      if (!additive_cv && nu_step_size <= 1)
      {
        writeLog_("Step size of nu <= 1 and additive_cv is false. Aborting!");
        return ILLEGAL_PARAMETERS;
      }
      double nu_stop = getDoubleOption_("cv:nu_stop");

      start_values.insert(make_pair(SVMWrapper::NU, nu_start));
      step_sizes.insert(make_pair(SVMWrapper::NU, nu_step_size));
      end_values.insert(make_pair(SVMWrapper::NU, nu_stop));
    }

    border_length = getIntOption_("border_length");
    svm.setParameter(SVMWrapper::BORDER_LENGTH, border_length);

    sigma = getDoubleOption_("sigma");
    svm.setParameter(SVMWrapper::SIGMA, sigma);

    k_mer_length = getIntOption_("k_mer_length");

    sigma_start = 0.;
    sigma_step_size = 0.;
    sigma_stop = 0.;
    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO
       && !skip_cv)
    {
      sigma_start = getDoubleOption_("cv:sigma_start");
      sigma_step_size = getDoubleOption_("cv:sigma_step_size");
      if (!additive_cv && sigma_step_size <= 1)
      {
        writeLog_("Step size of sigma <= 1 and additive_cv is false. Aborting!");
        return ILLEGAL_PARAMETERS;
      }
      sigma_stop = getDoubleOption_("cv:sigma_stop");

      start_values.insert(make_pair(SVMWrapper::SIGMA, sigma_start));
      step_sizes.insert(make_pair(SVMWrapper::SIGMA, sigma_step_size));
      end_values.insert(make_pair(SVMWrapper::SIGMA, sigma_stop));

      debug_string = "CV from sigma = " + String(sigma_start) +
                     " to sigma = " + String(sigma_stop) + " with step size " +
                     String(sigma_step_size);
      writeDebug_(debug_string, 1);
    }

    if (!skip_cv && !start_values.empty())
    {
      number_of_runs = getIntOption_("cv:number_of_runs");
      writeDebug_(String("Number of CV runs: ") + String(number_of_runs), 1);

      number_of_partitions = getIntOption_("cv:number_of_partitions");
      writeDebug_(String("Number of CV partitions: ") + String(number_of_partitions), 1);

      additive_cv = getFlag_("additive_cv");
    }

    Int debug_level = getIntOption_("debug");
    non_redundant = !(getFlag_("redundant"));

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    String document_id;
    IdXMLFile().load(inputfile_positives, protein_identifications, identifications, document_id);
    IdXMLFile().load(inputfile_negatives, protein_identifications_negative, identifications_negative, document_id);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    for (Size i = 0; i < identifications.size(); i++)
    {
      const vector<PeptideHit>& temp_peptide_hits = identifications[i].getHits();
      Size temp_size = temp_peptide_hits.size();
      if (temp_size > 0)
      {
        for (Size j = 0; j < temp_size; ++j)
        {
          temp_peptide_hit = temp_peptide_hits[j];
          temp_string = temp_peptide_hit.getSequence().toUnmodifiedString();
          if (!non_redundant
             || find(training_peptides.begin(), training_peptides.end(), temp_string) == training_peptides.end())
          {
            training_peptides.push_back(temp_peptide_hit.getSequence().toUnmodifiedString());
          }
        }
      }
    }
    training_labels.resize(training_peptides.size(), 1.0);
    debug_string = String(training_labels.size()) + " positive sequences read";
    writeDebug_(debug_string, 1);

    if (training_peptides.size() > max_positive_count)
    {
      random_shuffle(training_peptides.begin(), training_peptides.end());
      training_peptides.resize(max_positive_count, "");
      training_labels.resize(max_positive_count, 1.);
    }
    debug_string = String(training_peptides.size()) + " positive sequences for training";
    writeDebug_(debug_string, 1);

    UInt counter = 0;

    vector<String> temp_training_peptides;
    for (Size i = 0; i < identifications_negative.size(); i++)
    {
      const vector<PeptideHit>& temp_peptide_hits = identifications_negative[i].getHits();
      Size temp_size = temp_peptide_hits.size();
      if (temp_size > 0)
      {
        for (Size j = 0; j < temp_size; ++j)
        {
          temp_peptide_hit = temp_peptide_hits[j];
          temp_string = temp_peptide_hit.getSequence().toUnmodifiedString();
          if (find(training_peptides.begin(), training_peptides.end(), temp_string) != training_peptides.end())
          {
            writeLog_("Peptides are not allowed to occur in the positive and the negative set. Example: '" + temp_string + "'");
            return ILLEGAL_PARAMETERS;
          }

          if (!non_redundant
             || find(training_peptides.begin(), training_peptides.end(), temp_string) == training_peptides.end())
          {
            temp_training_peptides.push_back(temp_peptide_hit.getSequence().toUnmodifiedString());
            training_labels.push_back(-1.0);
            ++counter;
          }
        }
      }
    }
    if (non_redundant)
    {
      debug_string = String(counter) + " non redundant negative sequences read";
    }
    else
    {
      debug_string = String(counter) + " negative sequences read";
    }
    writeDebug_(debug_string, 1);
    if (temp_training_peptides.size() > max_negative_count)
    {
      random_shuffle(temp_training_peptides.begin(), temp_training_peptides.end());
      temp_training_peptides.resize(max_negative_count, "");
      training_labels.resize(training_peptides.size() + max_negative_count, -1.);
    }
    training_peptides.insert(training_peptides.end(),
                             temp_training_peptides.begin(),
                             temp_training_peptides.end());

    debug_string = String(temp_training_peptides.size()) + " negative sequences for training";
    writeDebug_(debug_string, 1);
    temp_training_peptides.clear();

    if (temp_type == LINEAR || temp_type == POLY || temp_type == RBF)
    {
      UInt maximum_sequence_length = 50;
      encoded_training_sample =
        encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(training_peptides,
                                                                   training_labels,
                                                                   allowed_amino_acid_characters,
                                                                   maximum_sequence_length);
    }
    else if (temp_type == SVMWrapper::OLIGO)
    {
      encoded_training_sample =
        encoder.encodeLibSVMProblemWithOligoBorderVectors(training_peptides,
                                                          training_labels,
                                                          k_mer_length,
                                                          allowed_amino_acid_characters,
                                                          svm.getIntParameter(SVMWrapper::BORDER_LENGTH));
    }

    if (!start_values.empty())
    {
      String digest = "";
      bool output_flag = false;
      if (debug_level >= 1)
      {
        output_flag = true;
        vector<String> parts;
        outputfile_name.split('/', parts);
        if (parts.empty())
        {
          digest = outputfile_name;
        }
        else
        {
          digest = parts[parts.size() - 1];
        }
      }
      SVMData dummy;
      double cv_quality = svm.performCrossValidation(encoded_training_sample,
                                                         dummy,
                                                         false,
                                                         start_values,
                                                         step_sizes,
                                                         end_values,
                                                         number_of_partitions,
                                                         number_of_runs,
                                                         optimized_parameters,
                                                         additive_cv,
                                                         output_flag,
                                                         "performances_" + digest + ".txt");

      String debug_string = "Best parameters found in cross validation:";

      for (parameters_iterator = optimized_parameters.begin();
           parameters_iterator != optimized_parameters.end();
           ++parameters_iterator)
      {
        svm.setParameter(parameters_iterator->first,
                         parameters_iterator->second);
        if (parameters_iterator->first == SVMWrapper::DEGREE)
        {
          debug_string += " degree: " + String(parameters_iterator->second);
        }
        else if (parameters_iterator->first == SVMWrapper::C)
        {
          debug_string += " C: " + String(parameters_iterator->second);
        }
        else if (parameters_iterator->first == SVMWrapper::NU)
        {
          debug_string += " nu: " + String(parameters_iterator->second);
        }
        else if (parameters_iterator->first == SVMWrapper::SIGMA)
        {
          debug_string += " sigma: " + String(parameters_iterator->second);
        }
      }
      debug_string += " with performance " + String(cv_quality);
      writeDebug_(debug_string, 1);
    }

    svm.train(encoded_training_sample);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    svm.saveModel(outputfile_name);

    // If the oligo-border kernel is used some additional information has to be stored.
    if (temp_type == SVMWrapper::OLIGO)
    {
      String outfile_name = getStringOption_("out");
      String param_outfile_name = getStringOption_("out_oligo_params");
      String trainset_outfile_name = getStringOption_("out_oligo_trainset");

      // Fallback to reasonable defaults if additional outfiles are not specified = empty.
      if (param_outfile_name.empty())
      {
        param_outfile_name = outfile_name + "_additional_parameters";
        writeLog_("Warning: Using OLIGO kernel but out_oligo_params was not specified. Trying to write to: " + param_outfile_name);
      }

      if (trainset_outfile_name.empty())
      {
        trainset_outfile_name = outfile_name + "_samples";
        writeLog_("Warning: Using OLIGO kernel but out_oligo_trainset was not specified. Trying to write to: " + trainset_outfile_name);
      }
      encoder.storeLibSVMProblem(trainset_outfile_name, encoded_training_sample);
      additional_parameters.setValue("kernel_type", temp_type);

      if (temp_type == SVMWrapper::OLIGO)
      {
        additional_parameters.setValue("border_length", svm.getIntParameter(SVMWrapper::BORDER_LENGTH));
        additional_parameters.setValue("k_mer_length", k_mer_length);
        additional_parameters.setValue("sigma", svm.getDoubleParameter(SVMWrapper::SIGMA));
      }

      ParamXMLFile paramFile;
      paramFile.store(param_outfile_name, additional_parameters);
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPPTModel tool;
  return tool.main(argc, argv);
}

/// @endcond
