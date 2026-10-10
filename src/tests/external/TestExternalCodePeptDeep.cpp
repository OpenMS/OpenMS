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
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#if defined(_MSC_VER) && defined(_DEBUG)
#include <crtdbg.h>
#endif

// Loads the models and checks the RT prediction; throws what the library throws.
static int run(const std::string& models)
{
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

int main(int argc, char** argv)
{
#if defined(_MSC_VER) && defined(_DEBUG)
  // The Debug C runtime reports a failed assertion and abort() through a dialog box
  // and waits for a click, so an uncaught exception or an assertion in the library
  // or the runtime would keep this program alive until the test timeout kills it.
  // Report to stderr instead and exit.
  _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
  _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);
#endif
  if (argc != 2)
  {
    std::cerr << "usage: TestExternalCodePeptDeep <directory with the PeptDeep models>\n";
    return 2;
  }
  try
  {
    return run(std::string(argv[1]) + "/");
  }
  catch (const std::exception& e)
  {
    std::cerr << "TestExternalCodePeptDeep: " << e.what() << "\n";
    return 1;
  }
}
