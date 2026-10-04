// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Consumer check of a WITH_ONNX installation: the PeptDeep inference classes,
// the models the installation ships and its ONNX Runtime work together. Every
// model is loaded, and the RT model has to reproduce the prediction of the
// AlphaPeptDeep reference (rt_pred_onnx of
// src/tests/class_tests/openms/data/peptdeep_irt_peptides_predicted.csv).

#include <OpenMS/ML/PEPTDEEP/PeptDeepCCSInference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "usage: TestExternalCodePeptDeep <directory with the PeptDeep models>\n";
    return 2;
  }
  const std::string models = std::string(argv[1]) + "/";

  OpenMS::PeptDeepCCSInference ccs(models + "peptdeep_ccs_dynamic.onnx");
  OpenMS::PeptDeepMS2Inference ms2(models + "peptdeep_ms2_dynamic.onnx");
  OpenMS::PeptDeepRTInference rt(models + "peptdeep_rt_dynamic.onnx");

  const std::vector<std::string> peptides = {"LGGNEQVTR", "GAGSSEPVTGLDAK"};
  const float expected[] = {0.07280362f, 0.2711957f};
  const std::vector<float> predicted = rt.predictRT(peptides);
  if (predicted.size() != peptides.size())
  {
    std::cerr << "predictRT returned " << predicted.size() << " values for " << peptides.size() << " peptides\n";
    return 1;
  }
  int rc = 0;
  for (size_t i = 0; i < predicted.size(); ++i)
  {
    const float diff = std::fabs(predicted[i] - expected[i]);
    std::cout << peptides[i] << ": " << predicted[i] << " (expected " << expected[i] << ")\n";
    if (diff > 1e-4f)
    {
      std::cerr << "RT prediction for " << peptides[i] << " differs from the reference by " << diff << "\n";
      rc = 1;
    }
  }
  return rc;
}
