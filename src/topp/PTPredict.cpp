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
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <map>
#include <iterator>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_PTPredict PTPredict

    @brief This application is used to predict the likelihood
                 of peptides to be proteotypic.

    This method has been described in the publication

    Ole Schulz-Trieglaff, Nico Pfeifer, Clemens Gr&ouml;pl, Oliver Kohlbacher and Knut Reinert LC-MSsim - a simulation software for Liquid ChromatographyMass Spectrometry data
  BMC Bioinformatics 2008, 9:423.

    The input of this application is an svm model and an idXML
    file with peptide identifications. The svm model file is specified
    by the <b>svm_model</b> parameter in the command line or the ini file.
    This file should have been produced by the @ref TOPP_PTModel application.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_PTPredict.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_PTPredict.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPTPredict :
  public TOPPBase
{
public:
  TOPPPTPredict() :
    TOPPBase("PTPredict", "predicts the likelihood of peptides to be proteotypic via svm_model which is trained by PTModel")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file\n");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerInputFile_("svm_model", "<file>", "", "svm model in libsvm format (can be produced by PTModel)");
    setValidFormats_("svm_model", ListUtils::create<String>("txt"));
    registerIntOption_("max_number_of_peptides", "<int>", 100000, "the maximum number of peptides considered at once (bigger number will lead to faster results but needs more memory).\n", false);
  }

  ExitCodes main_(int, const char**) override
  {
    IdXMLFile idXML_file;
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    vector<String> peptides;
    vector<PeptideHit> temp_peptide_hits;
    SVMWrapper svm;
    LibSVMEncoder encoder;
    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
    vector<double> predicted_likelihoods;
    vector<double> predicted_labels;
    map<String, double> predicted_data;
    svm_problem* training_data = nullptr;
    UInt border_length = 0;
    UInt k_mer_length = 0;
    double sigma = 0;
    String temp_string = "";
    UInt maximum_length = 50;
    UInt max_number_of_peptides = getIntOption_("max_number_of_peptides");


    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String inputfile_name = getStringOption_("in");
    String outputfile_name = getStringOption_("out");

    String svmfile_name = getStringOption_("svm_model");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    svm.loadModel(svmfile_name);

    // Since the POBK is not included in the libsvm we have to load
    // additional parameters from additional files.
    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
    {
      inputFileReadable_(svmfile_name + "_additional_parameters", "svm_model (derived)");

      Param additional_parameters;
      ParamXMLFile paramFile;
      paramFile.load(svmfile_name + "_additional_parameters", additional_parameters);
      if (additional_parameters.getValue("kernel_type") != DataValue::EMPTY)
      {
        svm.setParameter(SVMWrapper::KERNEL_TYPE, ((String) additional_parameters.getValue("kernel_type")).toInt());
      }

      if (additional_parameters.getValue("border_length") == DataValue::EMPTY
         && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        writeLog_("No border length saved in additional parameters file. Aborting!");
        cout << "No border length saved in additional parameters file. Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }
      border_length = ((String)additional_parameters.getValue("border_length")).toInt();
      if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
         && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        writeLog_("No k-mer length saved in additional parameters file. Aborting!");
        cout << "No k-mer length saved in additional parameters file. Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }
      k_mer_length = ((String)additional_parameters.getValue("k_mer_length")).toInt();
      if (additional_parameters.getValue("sigma") == DataValue::EMPTY
         && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        writeLog_("No sigma saved in additional parameters file. Aborting!");
        cout << "No sigma saved in additional parameters file. Aborting!" << endl;
        return ILLEGAL_PARAMETERS;
      }
      sigma = ((String)additional_parameters.getValue("sigma")).toFloat();

    }
    String document_id;
    idXML_file.load(inputfile_name, protein_identifications, identifications, document_id);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    for (Size i = 0; i < identifications.size(); i++)
    {
      temp_peptide_hits = identifications[i].getHits();
      for (Size j = 0; j < temp_peptide_hits.size(); j++)
      {
        peptides.push_back(temp_peptide_hits[j].getSequence().toUnmodifiedString());
      }
    }

    vector<double> labels;
    labels.resize(peptides.size(), 0);

    vector<String>::iterator it_from = peptides.begin();
    vector<String>::iterator it_to = peptides.begin();
    while (it_from != peptides.end())
    {
      vector<String> temp_peptides;
      vector<double> temp_labels;
      UInt i = 0;
      while (i <= max_number_of_peptides && it_to != peptides.end())
      {
        ++it_to;
        ++i;
      }

      temp_peptides.insert(temp_peptides.end(), it_from, it_to);
      temp_labels.resize(temp_peptides.size(), 0);

      svm_problem* prediction_data = nullptr;

      if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) != SVMWrapper::OLIGO)
      {
        prediction_data =
          encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(temp_peptides,
                                                                     temp_labels,
                                                                     allowed_amino_acid_characters,
                                                                     maximum_length);
      }
      else
      {
        prediction_data = encoder.encodeLibSVMProblemWithOligoBorderVectors(temp_peptides,
                                                                            temp_labels,
                                                                            k_mer_length,
                                                                            allowed_amino_acid_characters,
                                                                            border_length);
      }

      if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        inputFileReadable_((svmfile_name + "_samples").c_str(), "svm_model (derived)");

        training_data = encoder.loadLibSVMProblem(svmfile_name + "_samples");
        svm.setTrainingSample(training_data);

        svm.setParameter(SVMWrapper::BORDER_LENGTH, (Int) border_length);
        svm.setParameter(SVMWrapper::SIGMA, sigma);
      }
      svm.getSVCProbabilities(prediction_data, predicted_likelihoods, predicted_labels);

      for (Size i = 0; i < temp_peptides.size(); i++)
      {
        predicted_data.insert(make_pair(temp_peptides[i],
                                        (predicted_likelihoods[i])));
      }
      predicted_likelihoods.clear();
      predicted_labels.clear();
      LibSVMEncoder::destroyProblem(prediction_data);

      it_from = it_to;
    }

    for (Size i = 0; i < identifications.size(); i++)
    {
      temp_peptide_hits = identifications[i].getHits();
      for (Size j = 0; j < temp_peptide_hits.size(); j++)
      {
        double temp_likelihood = predicted_data[temp_peptide_hits[j].getSequence().toUnmodifiedString()];

        temp_peptide_hits[j].setMetaValue("predicted_PT", temp_likelihood);
      }
      identifications[i].setHits(temp_peptide_hits);
    }
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    idXML_file.store(outputfile_name,
                     protein_identifications,
                     identifications);
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPPTPredict tool;
  return tool.main(argc, argv);
}

/// @endcond
