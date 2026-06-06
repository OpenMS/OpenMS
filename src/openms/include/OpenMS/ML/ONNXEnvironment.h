#pragma once

#include <onnxruntime_cxx_api.h>

namespace OpenMS {

// Using a C++17 inline function ensures all classes share this exact same instance
inline Ort::Env& getONNXEnvironment() {
    static Ort::Env global_env{ORT_LOGGING_LEVEL_WARNING, "OpenMS_ONNX"};
    return global_env;
}

} // namespace OpenMS