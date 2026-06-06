#pragma once

#include <string>
#include <vector>
#include <memory>
#include <onnxruntime_cxx_api.h>
#include <OpenMS/config.h>

namespace OpenMS {

/// @brief Inference engine for PeptDeep MS2 Fragment Intensity prediction.
///
/// This class parses and tokenizes raw peptide string sequences into structural integer arrays
/// and evaluates them through the ONNX Runtime backend to predict fragment intensities.
class OPENMS_DLLAPI PeptDeepMS2Inference {
public:
    /// @brief Constructor initializes the ONNX session with the specified model graph.
    /// @param model_path Absolute or relative path to the PeptDeep MS2 ONNX model file.
    /// @throws Exception::BaseException if the model file cannot be found or loaded by the runtime.
    PeptDeepMS2Inference(const std::string& model_path);

    /// @brief Predicts MS2 fragment intensities for a single peptide sequence.
    /// @param peptide_sequence The raw uppercase peptide string (e.g., "PEPTIDEK").
    /// @param charge The precursor charge state.
    /// @param nce The normalized collision energy (NCE) applied during fragmentation.
    /// @param instrument_index Categorical integer representing the MS instrument type (default: 0).
    /// @return A flattened vector of predicted float fragment intensities representing the generated spectra.
    /// @throws std::invalid_argument if the sequence string is empty or contains unmapped amino acid characters.
    std::vector<float> predictMS2(const std::string& peptide_sequence,
                                  float charge,
                                  float nce,
                                  int64_t instrument_index = 0);

private:
    Ort::SessionOptions session_options_;
    std::unique_ptr<Ort::Session> session_;
    Ort::MemoryInfo memory_info_ = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
};

} // namespace OpenMS