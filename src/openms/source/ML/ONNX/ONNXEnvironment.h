// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#pragma once

// PRIVATE OpenMS-internal header. Intentionally NOT under include/OpenMS/ and
// NOT installed, so the ONNX Runtime third-party headers and the Ort::
// namespace never leak into the public OpenMS API. Only .cpp files in
// source/ML/ may include it.
#include <onnxruntime_cxx_api.h>

namespace OpenMS
{
    /**
   * @brief Returns the process-wide singleton ONNX Runtime environment.
   * * Thread-safe; initialized exactly once per application lifecycle.
   */
  Ort::Env& getONNXEnvironment();
}