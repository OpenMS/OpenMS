#pragma once

#include <onnxruntime_cxx_api.h>

namespace OpenMS
{
    // Forward declaration of the singleton environment.
    Ort::Env& getONNXEnvironment();
}