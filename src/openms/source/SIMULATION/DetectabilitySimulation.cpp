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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/DetectabilitySimulation.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

namespace OpenMS
{

  DetectabilitySimulation::DetectabilitySimulation() :
    DefaultParamHandler("DetectabilitySimulation")
  {
    setDefaultParams_();
  }

  DetectabilitySimulation::DetectabilitySimulation(const DetectabilitySimulation& source) :
    DefaultParamHandler(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  DetectabilitySimulation& DetectabilitySimulation::operator=(const DetectabilitySimulation& source)
  {
    setParameters(source.getParameters());
    updateMembers_();
    return *this;
  }

  DetectabilitySimulation::~DetectabilitySimulation()
  {
  }

  void DetectabilitySimulation::filterDetectability(SimTypes::FeatureMapSim& features)
  {
    LOG_INFO << "Detectability Simulation ... started" << std::endl;
    if (param_.getValue("dt_simulation_on") == "true")
    {
      svmFilter_(features);
    }
    else
    {
      noFilter_(features);
    }
  }

  void DetectabilitySimulation::noFilter_(SimTypes::FeatureMapSim& features)
  {
    // set detectibility to 1.0 for all given peptides
    double defaultDetectibility = 1.0;

    for (SimTypes::FeatureMapSim::iterator feature_it = features.begin();
         feature_it != features.end();
         ++feature_it)
    {
      (*feature_it).setMetaValue("detectability", defaultDetectibility);
    }
  }

  void DetectabilitySimulation::predictDetectabilities(vector<String>& peptides_vector, vector<double>& labels,
                                                       vector<double>& detectabilities)
  {
    // The support vector machine
    SVMWrapper svm;

    // initialize support vector machine
    LibSVMEncoder encoder;
    UInt k_mer_length = 0;
    double sigma = 0.0;
    UInt border_length = 0;

    if (File::readable(dt_model_file_))
    {
      svm.loadModel(dt_model_file_);
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "DetectibilitySimulation got invalid parameter. 'dt_model_file' " + dt_model_file_ + " is not readable");
    }

    // load additional parameters
    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
    {
      String add_paramfile = dt_model_file_ + "_additional_parameters";
      if (!File::readable(add_paramfile))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "DetectibilitySimulation: SVM parameter file " + add_paramfile + " is not readable");
      }

      Param additional_parameters;
      ParamXMLFile paramFile;
      paramFile.load(add_paramfile, additional_parameters);

      if (additional_parameters.getValue("border_length") == DataValue::EMPTY
         && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "DetectibilitySimulation: No border length defined in additional parameters file.");
      }
      border_length = ((String)additional_parameters.getValue("border_length")).toInt();
      if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
         && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "DetectibilitySimulation: No k-mer length defined in additional parameters file.");
      }
      k_mer_length = ((String)additional_parameters.getValue("k_mer_length")).toInt();

      if (additional_parameters.getValue("sigma") == DataValue::EMPTY
         && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "DetectibilitySimulation: No sigma defined in additional parameters file.");
      }

      sigma = ((String)additional_parameters.getValue("sigma")).toFloat();
    }

    if (File::readable(dt_model_file_))
    {
      svm.setParameter(SVMWrapper::BORDER_LENGTH, (Int) border_length);
      svm.setParameter(SVMWrapper::SIGMA, sigma);
      // to obtain probabilities
      svm.setParameter(SVMWrapper::PROBABILITY, 1);
    }
    // loading training data
    String sample_file = dt_model_file_ + "_samples";
    svm_problem* training_data = nullptr;
    if (File::readable(sample_file))
    {
      training_data = encoder.loadLibSVMProblem(sample_file);
      svm.setTrainingSample(training_data);
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "DetectibilitySimulation: SVM sample file " + sample_file + " is not readable");
    }


    LOG_INFO << "Predicting peptide detectabilities..    " << endl;

    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";

    // Encoding test data
    vector<double> probs;
    probs.resize(peptides_vector.size(), 0);

    svm_problem* prediction_data = encoder.encodeLibSVMProblemWithOligoBorderVectors(peptides_vector, probs,
                                                                                     k_mer_length,
                                                                                     allowed_amino_acid_characters,
                                                                                     svm.getIntParameter(SVMWrapper::BORDER_LENGTH));

    svm.getSVCProbabilities(prediction_data, detectabilities, labels);

    // clean up when finished with prediction
    delete prediction_data;
    delete training_data;
  }

  void DetectabilitySimulation::svmFilter_(SimTypes::FeatureMapSim& features)
  {

    // transform featuremap to peptides vector
    vector<String> peptides_vector(features.size());
    for (Size i = 0; i < features.size(); ++i)
    {
      peptides_vector[i] = features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString();
    }

    vector<double> labels;
    vector<double> detectabilities;
    predictDetectabilities(peptides_vector, labels, detectabilities);


    // copy all meta data stored in the feature map
    SimTypes::FeatureMapSim temp_copy(features);
    temp_copy.clear(false);

    for (Size i = 0; i < peptides_vector.size(); ++i)
    {

      if (detectabilities[i] > min_detect_)
      {
        features[i].setMetaValue("detectability", detectabilities[i]);
        temp_copy.push_back(features[i]);
      }
#ifdef DEBUG_SIM
      cout << detectabilities[i] << " " << min_detect_ << endl;
#endif
    }

    features.swap(temp_copy);
  }

  void DetectabilitySimulation::setDefaultParams_()
  {
    defaults_.setValue("dt_simulation_on", "false", "Modelling detectibility enabled? This can serve as a filter to remove peptides which ionize badly, thus reducing peptide count");
    defaults_.setValidStrings("dt_simulation_on", ListUtils::create<String>("true,false"));
    defaults_.setValue("min_detect", 0.5, "Minimum peptide detectability accepted. Peptides with a lower score will be removed");
    defaults_.setValue("dt_model_file", "examples/simulation/DTPredict.model", "SVM model for peptide detectability prediction");
    defaultsToParam_();
  }

  void DetectabilitySimulation::updateMembers_()
  {
    min_detect_ = param_.getValue("min_detect");
    dt_model_file_ = param_.getValue("dt_model_file");
    if (!File::readable(dt_model_file_)) // look in OPENMS_DATA_PATH
    {
      dt_model_file_ = File::find(dt_model_file_);
    }
  }

}
