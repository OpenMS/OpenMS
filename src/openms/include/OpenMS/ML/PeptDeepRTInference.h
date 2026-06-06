#pragma once

#include <string>
#include <vector>
#include <memory>

// Include the ONNX Runtime C++ API
#include <onnxruntime_cxx_api.h>

namespace OpenMS
{
    class __attribute__((visibility("default"))) PeptDeepRTInference
    {
    public:
        /**
         * @brief Constructor initializes the ONNX environment and loads the model
         * @param model_path Absolute path to peptdeep_rt_dynamic.onnx
         */
        explicit PeptDeepRTInference(const std::string& model_path);

        /**
         * @brief Destructor
         */
        ~PeptDeepRTInference();

        /**
         * @brief Predicts Retention Times for a list of peptide sequences
         * @param peptides A vector of peptide strings (e.g., "PEPTIDEK")
         * @return A vector of predicted RT values corresponding to the input peptides
         */
        std::vector<float> predictRT(const std::vector<std::string>& peptides);

    private:
        // ONNX Runtime components required to keep the model loaded in memory
        Ort::SessionOptions session_options_;
        std::unique_ptr<Ort::Session> session_;
        Ort::AllocatorWithDefaultOptions allocator_;

        // Helper function to convert text peptides into integer tensors for PeptDeep
        std::vector<int64_t> tokenizePeptides(const std::vector<std::string>& peptides, size_t& max_length);
    };
} // namespace OpenMS