#include <OpenMS/ML/PeptDeepMS2Inference.h>
#include <OpenMS/ML/ONNXEnvironment.h>
#include <onnxruntime_cxx_api.h>
#include <stdexcept>

namespace OpenMS {

namespace {
    inline int64_t getAAIndex(char aa) {
        if (aa >= 'A' && aa <= 'Z') return aa - 'A' + 1;
        if (aa >= 'a' && aa <= 'z') return aa - 'a' + 1;
        return 0; // 0 serves as the padding and unknown token
    }
}

// 1. THE HIDDEN IMPLEMENTATION STRUCT
struct PeptDeepMS2Inference::Impl
{
    Ort::SessionOptions session_options_;
    std::unique_ptr<Ort::Session> session_;
    Ort::MemoryInfo memory_info_ = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

    // Constructor
    Impl(const std::string& model_path)
    {
        try
        {
            session_options_.SetIntraOpNumThreads(1);
            session_options_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
            session_ = std::make_unique<Ort::Session>(getONNXEnvironment(), model_path.c_str(), session_options_);
        }
        catch (const Ort::Exception& e)
        {
            throw std::runtime_error(std::string("ONNX Runtime initialization failed: ") + e.what());
        }
    }

    // Core Prediction Logic
    std::vector<float> predictMS2(const std::string& peptide_sequence, float charge, float nce, int64_t instrument_index)
    {
        int64_t seq_length = peptide_sequence.size();
        int64_t batch_size = 1;

        if (seq_length == 0) {
            throw std::invalid_argument("Peptide sequence length cannot be zero.");
        }

        std::vector<int64_t> aa_indices;
        aa_indices.reserve(seq_length);

        for (char aa : peptide_sequence) {
            aa_indices.push_back(getAAIndex(aa));
        }

        std::vector<int64_t> aa_shape = {batch_size, seq_length};
        Ort::Value aa_tensor = Ort::Value::CreateTensor<int64_t>(
            memory_info_, const_cast<int64_t*>(aa_indices.data()), aa_indices.size(), aa_shape.data(), aa_shape.size());

        // Query mod_x dimensions dynamically instead of hardcoding 109
        Ort::TypeInfo mod_type_info = session_->GetInputTypeInfo(1);
        auto mod_tensor_info = mod_type_info.GetTensorTypeAndShapeInfo();
        std::vector<int64_t> mod_shape = mod_tensor_info.GetShape();

        if (mod_shape.size() < 3) {
            throw std::runtime_error("PeptDeep MS2 model input 'mod_x' must have 3 dimensions [batch, length, features].");
        }

        mod_shape[0] = batch_size;
        mod_shape[1] = seq_length;

        int64_t total_mod_elements = 1;
        for (int64_t dim : mod_shape) {
            if (dim < 0) {
                total_mod_elements *= 109; // Safeguard dynamic feature dimensions
            } else {
                total_mod_elements *= dim;
            }
        }

        std::vector<float> mod_x_data(total_mod_elements, 0.0f);
        Ort::Value mod_tensor = Ort::Value::CreateTensor<float>(
            memory_info_, mod_x_data.data(), mod_x_data.size(), mod_shape.data(), mod_shape.size());

        std::vector<float> charge_data = {charge};
        std::vector<int64_t> charge_shape = {batch_size, 1};
        Ort::Value charge_tensor = Ort::Value::CreateTensor<float>(
            memory_info_, charge_data.data(), charge_data.size(), charge_shape.data(), charge_shape.size());

        std::vector<float> nce_data = {nce};
        std::vector<int64_t> nce_shape = {batch_size, 1};
        Ort::Value nce_tensor = Ort::Value::CreateTensor<float>(
            memory_info_, nce_data.data(), nce_data.size(), nce_shape.data(), nce_shape.size());

        std::vector<int64_t> inst_data = {instrument_index};
        std::vector<int64_t> inst_shape = {batch_size};
        Ort::Value inst_tensor = Ort::Value::CreateTensor<int64_t>(
            memory_info_, inst_data.data(), inst_data.size(), inst_shape.data(), inst_shape.size());

        const char* input_names[] = {"aa_indices", "mod_x", "charges", "nce", "instrument_indices"};
        std::vector<Ort::Value> input_tensors;
        input_tensors.push_back(std::move(aa_tensor));
        input_tensors.push_back(std::move(mod_tensor));
        input_tensors.push_back(std::move(charge_tensor));
        input_tensors.push_back(std::move(nce_tensor));
        input_tensors.push_back(std::move(inst_tensor));

        const char* output_names[] = {"intensities"};

        auto output_tensors = session_->Run(
            Ort::RunOptions{nullptr}, input_names, input_tensors.data(), input_tensors.size(), output_names, 1
        );

        float* floatarr = output_tensors.front().GetTensorMutableData<float>();
        size_t output_count = output_tensors.front().GetTensorTypeAndShapeInfo().GetElementCount();

        return std::vector<float>(floatarr, floatarr + output_count);
    }
};

// 2. THE PUBLIC WRAPPER (Forwarding calls to Impl)

PeptDeepMS2Inference::PeptDeepMS2Inference(const std::string& model_path)
    : pimpl_(std::make_unique<Impl>(model_path))
{}

PeptDeepMS2Inference::~PeptDeepMS2Inference() = default;

std::vector<float> PeptDeepMS2Inference::predictMS2(const std::string& peptide_sequence, float charge, float nce, int64_t instrument_index)
{
    return pimpl_->predictMS2(peptide_sequence, charge, nce, instrument_index);
}

} // namespace OpenMS