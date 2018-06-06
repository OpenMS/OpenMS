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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <map>
#include <numeric>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_RTModel RTModel

    @brief Used to train a model for peptide retention time prediction or peptide separation prediction.

    For retention time prediction, a support vector machine is
    trained with peptide sequences and their measured retention
    times.
    For peptide separation prediction, two files have to be given:
    One file contains the positive examples (the peptides which
    are collected) and the other contains the negative examples
    (the flowthrough peptides).

    These methods and applications of this model are described
    in the following publications:

    Nico Pfeifer, Andreas Leinenbach, Christian G. Huber and Oliver Kohlbacher
    Statistical learning of peptide retention behavior in chromatographic separations: A new kernel-based approach for computational proteomics.
    BMC Bioinformatics 2007, 8:468

    Nico Pfeifer, Andreas Leinenbach, Christian G. Huber and Oliver Kohlbacher
    Improving Peptide Identification in Proteome Analysis by a Two-Dimensional Retention Time Filtering Approach
    J. Proteome Res. 2009, 8(8):4109-15

    There are a number of parameters which
    can be changed for the svm (specified in the ini file and command line):
    <ul>
        <li>
            svm_type: the type of the svm (can be NU_SVR or
            EPSILON_SVR for RT prediction and is C_SVC for separation
            prediction)
        </li>
        <li>
            kernel_type: the kernel function (e.g., POLY for the
                polynomial kernel, LINEAR for the linear kernel or RBF for the gaussian kernel); we recommend
                SVMWrapper::OLIGO for our paired oligo-border kernel (POBK)
        </li>
        <li>
            border_length: border length for the POBK
        </li>
        <li>
            k_mer_length: length of the signals considered in the
            POBK
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
            nu: the nu parameter for nu-SVR
        </li>
        <li>
            p: the epsilon parameter for epsilon-SVR
        </li>
    </ul>

    <br>

    The last five parameters (sigma, degree, c, nu and p)
    can be used in a cross validation (CV) to find the best parameters according to the
    training set. Therefore you have to specify the start value of a
    parameter, the step size in which the parameters should be increased
    and a final value for the particular parameter such that the tested
    parameter is never bigger than the given final value. If you want
    to perform a cross validation for example for the parameter c, enable CV (across all 5 parameters) and
    set @em skip_cv to <b>false</b> in the INI file. This can be easily done with using the INIFileEditor.


    Furthermore, you can specify the number of partitions for the CV with
    <b>number_of_partitions</b> in the ini file and the number of runs
    with <b>number_of_runs</b>.

    <br>
    Consequently you have two choices to use this application:

    <ol>
        <li>
            Set the parameters of the svm: The RTModel application will train
            the svm with the training data and store the svm model
        </li>
        <li>
            Give a range of parameters for which a CV should be performed:
            The RTModel application will perform a CV to find the best
            parameter combination in the given range and afterwards train
            the svm with the best parameters and the whole training data.
            Then the model is stored.
        </li>
    </ol>

    <br>
    The model can be used in @ref TOPP_RTPredict, to predict retention times
    for peptides or peptide separation depending on how you trained
    the model.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_RTModel.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_RTModel.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPRTModel :
  public TOPPBase
{
public:
  TOPPRTModel() :
    TOPPBase("RTModel", "Trains a model for the retention time prediction of peptides from a training set.")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "This is the name of the input file (RT prediction). It is assumed that the file type is idXML. Alternatively you can provide a .txt file having a sequence and the corresponding rt per line.\n", false);
    setValidFormats_("in", ListUtils::create<String>("idXML,txt"));
    registerInputFile_("in_positive", "<file>", "", "input file with positive examples (peptide separation prediction)\n", false);
    setValidFormats_("in_positive", ListUtils::create<String>("idXML"));
    registerInputFile_("in_negative", "<file>", "", "input file with negative examples (peptide separation prediction)\n", false);
    setValidFormats_("in_negative", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file: the model in libsvm format");
    setValidFormats_("out", ListUtils::create<String>("txt"));
    registerOutputFile_("out_oligo_params", "<file>", "", "output file with additional model parameters when using the OLIGO kernel", false);
    setValidFormats_("out_oligo_params", ListUtils::create<String>("paramXML"));
    registerOutputFile_("out_oligo_trainset", "<file>", "", "output file with the used training dataset when using the OLIGO kernel", false);
    setValidFormats_("out_oligo_trainset", ListUtils::create<String>("txt"));
    registerStringOption_("svm_type", "<type>", "NU_SVR", "the type of the svm (NU_SVR or EPSILON_SVR for RT prediction, automatically set\nto C_SVC for separation prediction)\n", false);
    setValidStrings_("svm_type", ListUtils::create<String>("NU_SVR,NU_SVC,EPSILON_SVR,C_SVC"));
    registerDoubleOption_("nu", "<float>", 0.5, "the nu parameter [0..1] of the svm (for nu-SVR)", false);
    setMinFloat_("nu", 0);
    setMaxFloat_("nu", 1);
    registerDoubleOption_("p", "<float>", 0.1, "the epsilon parameter of the svm (for epsilon-SVR)", false);
    registerDoubleOption_("c", "<float>", 1, "the penalty parameter of the svm", false);
    registerStringOption_("kernel_type", "<type>", "OLIGO", "the kernel type of the svm", false);
    setValidStrings_("kernel_type", ListUtils::create<String>("LINEAR,RBF,POLY,OLIGO"));
    registerIntOption_("degree", "<int>", 1, "the degree parameter of the kernel function of the svm (POLY kernel)\n", false);
    setMinInt_("degree", 1);
    registerIntOption_("border_length", "<int>", 22, "length of the POBK", false);
    setMinInt_("border_length", 1);
    registerDoubleOption_("max_std", "<float>", 10, "max standard deviation for a peptide to be included (if there are several ones for one peptide string)(median is taken)", false);
    setMinFloat_("max_std", 0.);
    registerIntOption_("k_mer_length", "<int>", 1, "k_mer length of the POBK", false);
    setMinInt_("k_mer_length", 1);
    registerDoubleOption_("sigma", "<float>", 5, "sigma of the POBK", false);
    registerDoubleOption_("total_gradient_time", "<time>", 1, "the time (in seconds) of the gradient (only for RT prediction)", false);
    setMinFloat_("total_gradient_time", 0.00001);
    registerFlag_("first_dim_rt", "if set the model will be built for first_dim_rt");
    registerFlag_("additive_cv", "if the step sizes should be interpreted additively (otherwise the actual value is multiplied\nwith the step size to get the new value");

    addEmptyLine_();
    registerTOPPSubsection_("cv", "Parameters for the grid search / cross validation:");
    registerFlag_("cv:skip_cv", "Set to enable Cross-Validation or set to true if the model should just be trained with 1 set of specified parameters.");
    registerIntOption_("cv:number_of_runs", "<int>", 1, "number of runs for the CV (each run creates a new random partition of the data)", false);
    setMinInt_("cv:number_of_runs", 1);
    registerIntOption_("cv:number_of_partitions", "<int>", 10, "number of CV partitions", false);
    setMinInt_("cv:number_of_partitions", 2);

    registerIntOption_("cv:degree_start", "<int>", 1, "starting point of degree", false);
    setMinInt_("cv:degree_start", 1);
    registerIntOption_("cv:degree_step_size", "<int>", 2, "step size point of degree", false);
    registerIntOption_("cv:degree_stop", "<int>", 4, "stopping point of degree", false);

    registerDoubleOption_("cv:p_start", "<float>", 1, "starting point of p", false);
    registerDoubleOption_("cv:p_step_size", "<float>", 10, "step size point of p", false);
    registerDoubleOption_("cv:p_stop", "<float>", 1000, "stopping point of p", false);

    registerDoubleOption_("cv:c_start", "<float>", 1, "starting point of c", false);
    registerDoubleOption_("cv:c_step_size", "<float>", 10, "step size of c", false);
    registerDoubleOption_("cv:c_stop", "<float>", 1000, "stopping point of c", false);

    registerDoubleOption_("cv:nu_start", "<float>", 0.3, "starting point of nu", false);
    setMinFloat_("cv:nu_start", 0);
    setMaxFloat_("cv:nu_start", 1);
    registerDoubleOption_("cv:nu_step_size", "<float>", 1.2, "step size of nu", false);
    registerDoubleOption_("cv:nu_stop", "<float>", 0.7, "stopping point of nu", false);
    setMinFloat_("cv:nu_stop", 0);
    setMaxFloat_("cv:nu_stop", 1);

    registerDoubleOption_("cv:sigma_start", "<float>", 1, "starting point of sigma", false);
    registerDoubleOption_("cv:sigma_step_size", "<float>", 1.3, "step size of sigma", false);
    registerDoubleOption_("cv:sigma_stop", "<float>", 15, "stopping point of sigma", false);

  }

  void loadStringLabelLines_(String                     filename,
                             std::vector<String>& sequences,
                             std::vector<double>& labels)
  {
    TextFile text_file(filename.c_str(), true);
    std::vector<String> parts;
    labels.clear();

    TextFile::ConstIterator it = text_file.begin();
    while (it != text_file.end())
    {
      it->split(' ', parts);
      if (parts.size() == 2)
      {
        sequences.push_back(parts[0].trim());
        labels.push_back(parts[1].trim().toDouble());
        ++it;
      }
      else
      {
        it->split('\v', parts);
        if (parts.size() == 2)
        {
          sequences.push_back(parts[0].trim());
          labels.push_back(parts[1].trim().toDouble());
          ++it;
        }
        else
        {
          it->split('\t', parts);
          if (parts.size() == 2)
          {
            sequences.push_back(parts[0].trim());
            labels.push_back(parts[1].trim().toDouble());
            ++it;
          }
          else
          {
            String debug_string = "found line '" + *it + "' in file which is not of the form <string> <label>\n";
            writeDebug_(debug_string, 1);
            ++it;
          }
        }
      }
    }
  }

  ExitCodes main_(Int, const char**) override
  {
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    vector<ProteinIdentification> protein_identifications_negative;
    vector<PeptideIdentification> identifications_negative;
    vector<String> training_peptides;
    vector<AASequence> training_modified_peptides;
    vector<double> training_retention_times;
    PeptideHit temp_peptide_hit;
    SVMWrapper svm;
    svm.setLogType(log_type_);
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
    pair<double, double> sigmas;
    Int temp_type = POLY;
    String debug_string = "";
    double sigma = 0.1;
    UInt k_mer_length = 1;
    Int border_length = 0;
    bool separation_prediction = false;
    map<String, double> redundant_peptides;
    map<AASequence, double> redundant_modified_peptides;
    double max_std = 0.;
    bool textfile_input = false;
    SVMData training_sample;
    bool first_dim_rt = false;
    bool skip_cv = false;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_positives = getStringOption_("in_positive");
    String inputfile_negatives = "";
    String inputfile_name = "";
    if (inputfile_positives != "")
    {
      inputfile_negatives = getStringOption_("in_negative");
      if (inputfile_negatives != "")
      {
        separation_prediction = true;
      }
      else
      {
        writeLog_("Positive peptides for separation prediction set but no negative peptides. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
    }
    else
    {
      inputfile_name = getStringOption_("in");
      textfile_input = (FileHandler::getTypeByFileName(inputfile_name) == FileTypes::TXT);
    }
    String outputfile_name = getStringOption_("out");
    additive_cv = getFlag_("additive_cv");
    skip_cv = getFlag_("cv:skip_cv");
    if (skip_cv) LOG_INFO << "Cross-validation disabled!\n";
    else LOG_INFO << "Cross-validation enabled!\n";

    float total_gradient_time = getDoubleOption_("total_gradient_time");
    max_std = getDoubleOption_("max_std");
    if (!separation_prediction && total_gradient_time < 0)
    {
      writeLog_("No total gradient time given for RT prediction. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    //SVM type
    String type = getStringOption_("svm_type");
    if (type == "NU_SVR" && !separation_prediction)
    {
      svm.setParameter(SVMWrapper::SVM_TYPE, NU_SVR);
    }
    else if (type == "EPSILON_SVR" && !separation_prediction)
    {
      svm.setParameter(SVMWrapper::SVM_TYPE, EPSILON_SVR);
    }
    else if ((separation_prediction && type == "C_SVC")
            || separation_prediction)
    {
      svm.setParameter(SVMWrapper::SVM_TYPE, C_SVC);
    }
    else
    {
      writeLog_("Illegal SVM type given. SVM type has to be either "
                + String("NU_SVR or EPSILON_SVR for RT prediction and ")
                + "C_SVC for separation prediction. Aborting!");
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
    if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVR || svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVC)
    {
      svm.setParameter(SVMWrapper::NU, getDoubleOption_("nu"));
    }
    else if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == EPSILON_SVR)
    {
      svm.setParameter(SVMWrapper::P, getDoubleOption_("p"));
    }

    //grid search parameters
    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == POLY)
    {
      svm.setParameter(SVMWrapper::DEGREE, getIntOption_("degree"));

      if (!skip_cv)
      {
        UInt degree_start = getIntOption_("cv:degree_start");
        UInt degree_step_size = getIntOption_("cv:degree_step_size");
        if (!additive_cv && degree_step_size <= 1)
        {
          writeLog_("Step size of degree <= 1 and additive_cv is false. Aborting!");
          return ILLEGAL_PARAMETERS;
        }
        UInt degree_stop = getIntOption_("cv:degree_stop");

        start_values.insert(make_pair(SVMWrapper::DEGREE, degree_start));
        step_sizes.insert(make_pair(SVMWrapper::DEGREE, degree_step_size));
        end_values.insert(make_pair(SVMWrapper::DEGREE, degree_stop));
      }
    }

    if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == EPSILON_SVR && !skip_cv)
    {
      double p_start = getDoubleOption_("cv:p_start");
      double p_step_size = getDoubleOption_("cv:p_step_size");
      if (!additive_cv && p_step_size <= 1)
      {
        writeLog_("Step size of p <= 1 and additive_cv is false. Aborting!");
        return ILLEGAL_PARAMETERS;
      }
      double p_stop = getDoubleOption_("cv:p_stop");

      start_values.insert(make_pair(SVMWrapper::P, p_start));
      step_sizes.insert(make_pair(SVMWrapper::P, p_step_size));
      end_values.insert(make_pair(SVMWrapper::P, p_stop));
    }

    if (!skip_cv)
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

    if ((svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVR || svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVC)
       && !skip_cv)
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
    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
    {
      border_length = getIntOption_("border_length");
    }

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
    if (!start_values.empty())
    {
      number_of_runs = getIntOption_("cv:number_of_runs");
      writeDebug_(String("Number of CV runs: ") + String(number_of_runs), 1);

      number_of_partitions = getIntOption_("cv:number_of_partitions");
      writeDebug_(String("Number of CV partitions: ") + String(number_of_partitions), 1);
    }

    first_dim_rt = getFlag_("first_dim_rt");

    Int debug_level = getIntOption_("debug");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    if (!separation_prediction)
    {
      if (textfile_input)
      {
        loadStringLabelLines_(inputfile_name, training_peptides, training_retention_times);
        for (Size i = 0; i < training_peptides.size(); ++i)
        {
          if (temp_type == SVMWrapper::OLIGO)
          {
            redundant_modified_peptides.insert(make_pair(AASequence::fromString(training_peptides[i]),
                                                         training_retention_times[i]));
          }
          else
          {
            redundant_peptides.insert(make_pair(training_peptides[i], training_retention_times[i]));
          }
        }
        training_peptides.clear();
        training_retention_times.clear();
      }
      else
      {
        String document_id;
        IdXMLFile().load(inputfile_name, protein_identifications, identifications, document_id);
      }
    }
    else
    {
      String document_id;
      IdXMLFile().load(inputfile_positives, protein_identifications, identifications, document_id);
      IdXMLFile().load(inputfile_negatives, protein_identifications_negative, identifications_negative, document_id);
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    if (!textfile_input)
    {
      for (Size i = 0; i < identifications.size(); i++)
      {
        Size temp_size = identifications[i].getHits().size();
        if (temp_size > 0)
        {
          if (temp_size == 1)
          {
            temp_peptide_hit = identifications[i].getHits()[0];
            if (separation_prediction)
            {
              training_retention_times.push_back(1.0);
              if (temp_type == SVMWrapper::OLIGO)
              {
                training_modified_peptides.push_back(temp_peptide_hit.getSequence());
              }
              else
              {
                training_peptides.push_back(temp_peptide_hit.getSequence().toUnmodifiedString());
              }
            }
            else
            {
              if (first_dim_rt)
              {
                if (temp_type != SVMWrapper::OLIGO)
                {
                  redundant_peptides.insert(make_pair(temp_peptide_hit.getSequence().toUnmodifiedString(),
                                                      (double)(identifications[i].getMetaValue("first_dim_rt"))));
                }
                else
                {
                  redundant_modified_peptides.insert(make_pair(temp_peptide_hit.getSequence(),
                                                               (double)(identifications[i].getMetaValue("first_dim_rt"))));
                }
              }
              else
              {
                if (temp_type != SVMWrapper::OLIGO)
                {
                  redundant_peptides.insert(make_pair(temp_peptide_hit.getSequence().toUnmodifiedString(),
                                                      identifications[i].getRT()));
                }
                else
                {
                  redundant_modified_peptides.insert(make_pair(temp_peptide_hit.getSequence(),
                                                               identifications[i].getRT()));
                }
              }
            }
          }
          else
          {
            writeLog_("For one spectrum there should not be more than one peptide."
                      "Please use the IDFilter with the -best:strict option to achieve this. Aborting!");
            writeLog_("Hits: ");
            for (vector<PeptideHit>::const_iterator it = identifications[i].getHits().begin();
                 it != identifications[i].getHits().end();
                 ++it)
            {
              writeLog_(String(it->getSequence().toUnmodifiedString()) + " score: " + String(it->getScore()));
            }
            return INPUT_FILE_CORRUPT;
          }
        }
      }
    } // end ! textfile input

    // Getting a non redundant training set. If there are several copies of one peptide,
    // the standard deviation is calculated. If this std is less or equal to the
    // maximal allowed std 'max_std', the peptide is added to the training set
    // with the median as retention time. Unique peptides are added immediately.
    if (!separation_prediction && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
    {
      map<AASequence, double>::iterator it = redundant_modified_peptides.begin();

      double temp_median = 0;
      double temp_mean = 0;
      vector<double> temp_values;
      pair<map<AASequence, double>::iterator, map<AASequence, double>::iterator> it_pair;

      while (it != redundant_modified_peptides.end())
      {
        temp_values.clear();
        double temp_variance = 0;

        it_pair = redundant_modified_peptides.equal_range(it->first);
        for (it = it_pair.first; it != it_pair.second; ++it)
        {
          temp_values.push_back(it->second);
        }
        if (temp_values.size() == 1)
        {
          temp_median = temp_values[0];
          temp_mean = temp_values[0];
        }
        else
        {
          sort(temp_values.begin(), temp_values.end());
          if (temp_values.size() % 2 == 1)
          {
            temp_median = temp_values[temp_values.size() / 2];
          }
          else
          {
            temp_median = ((double) temp_values[temp_values.size() / 2]
                           + temp_values[temp_values.size() / 2 - 1]) / 2;
          }

          temp_mean = accumulate(temp_values.begin(), temp_values.end(), 0.) / temp_values.size();

          for (Size j = 0; j < temp_values.size(); ++j)
          {
            temp_variance += (temp_values[j] - temp_mean) * (temp_values[j] - temp_mean);
          }
          temp_variance /= temp_values.size();
        }
        if (sqrt(temp_variance) <= max_std)
        {
          training_modified_peptides.push_back(it_pair.first->first);
          training_retention_times.push_back(temp_median);
        }
        else
        {
          debug_string = "Did not take peptide " + it_pair.first->first.toString() + " for training because"
                         + " there were several copies and std was " + String(temp_median)
                         + " while " + String(max_std) + " was allowed.";
          writeDebug_(debug_string, 1);
        }
      }
    }

    if (!separation_prediction && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) != SVMWrapper::OLIGO)
    {
      map<String, double>::iterator it = redundant_peptides.begin();
      double temp_median = 0;
      double temp_mean = 0;
      vector<double> temp_values;
      pair<map<String, double>::iterator, map<String, double>::iterator> it_pair;

      while (it != redundant_peptides.end())
      {
        temp_values.clear();
        double temp_variance = 0;

        it_pair = redundant_peptides.equal_range(it->first);
        for (it = it_pair.first; it != it_pair.second; ++it)
        {
          temp_values.push_back(it->second);
        }
        if (temp_values.size() == 1)
        {
          temp_median = temp_values[0];
          temp_mean = temp_values[0];
        }
        else
        {
          sort(temp_values.begin(), temp_values.end());
          if (temp_values.size() % 2 == 1)
          {
            temp_median = temp_values[temp_values.size() / 2];
          }
          else
          {
            temp_median = ((double) temp_values[temp_values.size() / 2]
                           + temp_values[temp_values.size() / 2 - 1]) / 2;
          }

          temp_mean = accumulate(temp_values.begin(), temp_values.end(), 0.) / temp_values.size();

          for (Size j = 0; j < temp_values.size(); ++j)
          {
            temp_variance += (temp_values[j] - temp_mean) * (temp_values[j] - temp_mean);
          }
          temp_variance /= temp_values.size();
        }
        if (sqrt(temp_variance) <= max_std)
        {
          training_peptides.push_back(it_pair.first->first);
          training_retention_times.push_back(temp_median);
        }
        else
        {
          debug_string = "Did not take peptide " + it_pair.first->first + " for training because"
                         + " there were several copies and std was " + String(temp_median)
                         + " while " + String(max_std) + " was allowed.";
          writeDebug_(debug_string, 1);
        }
      }
    }

    // For separation prediction there are two files needed
    if (separation_prediction)
    {
      for (Size i = 0; i < identifications_negative.size(); i++)
      {
        Size temp_size = identifications_negative[i].getHits().size();
        if (temp_size > 0)
        {
          if (temp_size == 1)
          {
            temp_peptide_hit = identifications_negative[i].getHits()[0];
            if (temp_type == SVMWrapper::OLIGO)
            {
              training_modified_peptides.push_back(temp_peptide_hit.getSequence());
            }
            else
            {
              training_peptides.push_back(temp_peptide_hit.getSequence().toUnmodifiedString());
            }

            training_retention_times.push_back(-1.0);
          }
          else
          {
            writeLog_("For one spectrum there should not be more than one peptide."
                      "Please use the IDFilter with the -best:strict option to achieve this. Aborting!");
            writeLog_("Hits: ");
            for (vector<PeptideHit>::const_iterator it = identifications_negative[i].getHits().begin();
                 it != identifications_negative[i].getHits().end();
                 ++it)
            {
              writeLog_(String(it->getSequence().toUnmodifiedString()) + " score: " + String(it->getScore()));
            }
            return INPUT_FILE_CORRUPT;
          }
        }
      }
    }

    if (!separation_prediction)
    {
      for (Size i = 0; i < training_retention_times.size(); i++)
      {
        training_retention_times[i] = training_retention_times[i] / total_gradient_time;
      }
    }

    if (temp_type == LINEAR || temp_type == POLY || temp_type == RBF)
    {
      // TODO What happens if the sequence exceeds this size? No error, but does it impact performance?
      // Why this magic number?
      UInt maximum_sequence_length = 50;
      encoded_training_sample =
        encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(training_peptides,
                                                                   training_retention_times,
                                                                   allowed_amino_acid_characters,
                                                                   maximum_sequence_length);
    }
    else if (temp_type == SVMWrapper::OLIGO)
    {
      encoder.encodeProblemWithOligoBorderVectors(training_modified_peptides,
                                                  k_mer_length,
                                                  allowed_amino_acid_characters,
                                                  svm.getIntParameter(SVMWrapper::BORDER_LENGTH),
                                                  training_sample.sequences);
      training_sample.labels = training_retention_times;
    }

    if (!skip_cv && !start_values.empty())
    {
      String digest = "";
      bool output_flag = false;
      if (debug_level >= 1)
      {
        output_flag = true;
        vector<String> parts;
        inputfile_name.split('/', parts);
        if (parts.empty())
        {
          digest = inputfile_name;
        }
        else
        {
          digest = parts[parts.size() - 1];
        }
      }
      double cv_quality = 0.0;

      if (temp_type == SVMWrapper::OLIGO)
      {
        debug_string = String(training_sample.sequences.size()) + " sequences for training, "
                       + training_sample.labels.size() + " labels for training";
        writeDebug_(debug_string, 1);

        cv_quality = svm.performCrossValidation(encoded_training_sample,
                                                training_sample,
                                                true,
                                                start_values,
                                                step_sizes,
                                                end_values,
                                                number_of_partitions,
                                                number_of_runs,
                                                optimized_parameters,
                                                additive_cv,
                                                output_flag,
                                                "performances_" + digest + ".txt");
      }
      else
      {
        cv_quality = svm.performCrossValidation(encoded_training_sample,
                                                training_sample,
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
      }
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
        else if (parameters_iterator->first == SVMWrapper::P)
        {
          debug_string += " P: " + String(parameters_iterator->second);
        }
        else if (parameters_iterator->first == SVMWrapper::SIGMA)
        {
          debug_string += " sigma: " + String(parameters_iterator->second);
        }
      }
      debug_string += " with performance " + String(cv_quality);
      writeDebug_(debug_string, 1);
    }

    if (temp_type == SVMWrapper::OLIGO)
    {
      svm.train(training_sample);
    }
    else
    {
      svm.train(encoded_training_sample);
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    svm.saveModel(outputfile_name);

    // If the oligo-border kernel is used some additional information has to be stored
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
      training_sample.store(trainset_outfile_name);
      additional_parameters.setValue("kernel_type", temp_type);

      if (!separation_prediction)
      {
        svm.getSignificanceBorders(training_sample, sigmas);

        additional_parameters.setValue("sigma_0", sigmas.first);
        additional_parameters.setValue("sigma_max", sigmas.second);
        if (first_dim_rt)
        {
          additional_parameters.setValue("first_dim_rt", "true");
        }
      }
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
  TOPPRTModel tool;
  return tool.main(argc, argv);
}

/// @endcond
