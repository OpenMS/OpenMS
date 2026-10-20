// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include "ONNXEnvironment.h"

#if defined(__EMSCRIPTEN_PTHREADS__)
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

namespace OpenMS
{
  Ort::Env& getONNXEnvironment()
  {
    // A thread-safe static instance ensures the ONNX Runtime environment
    // is initialized exactly once per application lifecycle.
#if defined(__EMSCRIPTEN_PTHREADS__)
    // WebAssembly with threads: ONNX Runtime's sessions use the environment's global thread pools
    // there by default (as onnxruntime-web, which creates them the same way), so the environment
    // has to have them. One pool of as many threads as OpenMP computes with, shared by all
    // sessions; its threads are started with the environment.
    static Ort::Env env = []
    {
      Ort::ThreadingOptions threading;
#ifdef _OPENMP
      threading.SetGlobalIntraOpNumThreads(std::max(1, omp_get_max_threads()));
#else
      threading.SetGlobalIntraOpNumThreads(1);
#endif
      threading.SetGlobalInterOpNumThreads(1);
      return Ort::Env(threading, ORT_LOGGING_LEVEL_WARNING, "OpenMS_ONNX_Instance");
    }();
#else
    static Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "OpenMS_ONNX_Instance");
#endif
    return env;
  }
}
