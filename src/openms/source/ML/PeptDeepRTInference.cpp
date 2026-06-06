#include <OpenMS/ML/PeptDeepRTInference.h>
#include <OpenMS/ML/ONNXEnvironment.h>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

namespace OpenMS
{
    // 1. ENGINE STARTUP (Constructor)
    PeptDeepRTInference::PeptDeepRTInference(const std::string& model_path)
{
    session_options_.SetIntraOpNumThreads(1);
    session_options_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

    // Uses the global shared environment:
    session_ = std::make_unique<Ort::Session>(getONNXEnvironment(), model_path.c_str(), session_options_);
}

    PeptDeepRTInference::~PeptDeepRTInference() = default;

    // 2. DATA TRANSLATOR (Tokenization)
    std::vector<int64_t> PeptDeepRTInference::tokenizePeptides(const std::vector<std::string>& peptides, size_t& max_length)
    {
        std::unordered_map<char, int64_t> vocab = {
            {'A', 1}, {'R', 2}, {'N', 3}, {'D', 4}, {'C', 5},
            {'E', 6}, {'Q', 7}, {'G', 8}, {'H', 9}, {'I', 10},
            {'L', 11}, {'K', 12}, {'M', 13}, {'F', 14}, {'P', 15},
            {'S', 16}, {'T', 17}, {'W', 18}, {'Y', 19}, {'V', 20}
        };

        // FORCE the sequence length to exactly 132 to match PeptDeep's strict ONNX graph
        max_length = 132;

        std::vector<int64_t> flat_tokens;
        flat_tokens.reserve(peptides.size() * max_length);

        // Convert strings to integers and pad with 0s up to 132
        for (const auto& p : peptides) {
            for (size_t i = 0; i < max_length; ++i) {
                if (i < p.length() && vocab.count(p[i])) {
                    flat_tokens.push_back(vocab[p[i]]);
                } else {
                    flat_tokens.push_back(0); // Zero-Padding to reach 132
                }
            }
        }
        return flat_tokens;
    }

    // 3. THE PREDICTION LOGIC
    std::vector<float> PeptDeepRTInference::predictRT(const std::vector<std::string>& peptides)
    {
        if (peptides.empty()) return {};

        // 1. Tokenize the input strings
        size_t max_seq_len = 132;
        std::vector<int64_t> input_tokens = tokenizePeptides(peptides, max_seq_len);
        std::vector<int64_t> seq_shape = { static_cast<int64_t>(peptides.size()), static_cast<int64_t>(max_seq_len) };

        // 2. DYNAMIC 3D TENSOR: Ask the model what shape "mod_x" needs to be
        Ort::TypeInfo mod_type_info = session_->GetInputTypeInfo(1);
        auto mod_tensor_info = mod_type_info.GetTensorTypeAndShapeInfo();
        std::vector<int64_t> mod_shape = mod_tensor_info.GetShape();

        // The model usually returns [-1, -1, FEATURE_SIZE].
        // We replace the flexible '-1' dimensions with our actual Batch Size and Sequence Length.
        mod_shape[0] = peptides.size();
        mod_shape[1] = max_seq_len;

        // Calculate the total number of zeroes needed to fill this 3D cube
        int64_t total_mod_elements = 1;
        for (int64_t dim : mod_shape) {
            total_mod_elements *= dim;
        }

        // Create the 3D zero-filled modifications array (0.0 means no modifications)
        std::vector<float> mod_x_data(total_mod_elements, 0.0f);

        // 3. Create the ONNX Memory Allocator
        Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

        // 4. Push BOTH tensors into a C++ Vector
        std::vector<Ort::Value> input_tensors;

        // Tensor 1: The Sequence (2D)
        input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
            memory_info, input_tokens.data(), input_tokens.size(), seq_shape.data(), seq_shape.size()
        ));

        // Tensor 2: The Modifications (3D)
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info, mod_x_data.data(), mod_x_data.size(), mod_shape.data(), mod_shape.size()
        ));

        // 5. Ask the model for the actual names of Door 1 (Sequence) and Door 2 (Mods)
        Ort::AllocatorWithDefaultOptions ort_alloc;
        Ort::AllocatedStringPtr seq_name_ptr = session_->GetInputNameAllocated(0, ort_alloc);
        Ort::AllocatedStringPtr mod_name_ptr = session_->GetInputNameAllocated(1, ort_alloc);
        Ort::AllocatedStringPtr output_name_ptr = session_->GetOutputNameAllocated(0, ort_alloc);

        const char* input_names[] = { seq_name_ptr.get(), mod_name_ptr.get() };
        const char* output_names[] = { output_name_ptr.get() };

        // 6. Fire the Engine!
        auto output_tensors = session_->Run(
            Ort::RunOptions{nullptr},
            input_names,
            input_tensors.data(), 2,
            output_names, 1
        );

        // 7. Extract the raw float predictions and return them
        float* floatarr = output_tensors.front().GetTensorMutableData<float>();
        std::vector<float> results(floatarr, floatarr + peptides.size());

        return results;
    }

} // namespace OpenMS