// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/ONNXEnvironment.h>

namespace OpenMS
{
  Ort::Env& getONNXEnvironment()
  {
    // A thread-safe static instance ensures the ONNX Runtime environment
    // is initialized exactly once per application lifecycle.
    static Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "OpenMS_ONNX_Instance");
    return env;
  }
}