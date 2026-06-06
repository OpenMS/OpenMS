#include <OpenMS/ML/PeptDeepMS2Inference.h>
#include <OpenMS/ML/ONNXEnvironment.h>
#include <stdexcept>
#include <map>

namespace OpenMS {

// The official PeptDeep 1-based token alphabet mapping
static const std::map<char, int64_t> aa_vocab = {
    {'A', 1}, {'R', 2}, {'N', 3}, {'D', 4}, {'C', 5},
    {'Q', 6}, {'E', 7}, {'G', 8}, {'H', 9}, {'I', 10},
    {'L', 11}, {'K', 12}, {'M', 13}, {'F', 14}, {'P', 15},
    {'S', 16}, {'T', 17}, {'W', 18}, {'Y', 19}, {'V', 20}
};

PeptDeepMS2Inference::PeptDeepMS2Inference(const std::string& model_path) {
    session_options_.SetIntraOpNumThreads(1);
    session_options_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    session_ = std::make_unique<Ort::Session>(getONNXEnvironment(), model_path.c_str(), session_options_);
}

std::vector<float> PeptDeepMS2Inference::predictMS2(const std::string& peptide_sequence,
                                                    float charge,
                                                    float nce,
                                                    int64_t instrument_index) {

    int64_t seq_length = peptide_sequence.size();
    int64_t batch_size = 1;

    if (seq_length == 0) {
        throw std::invalid_argument("Peptide sequence length cannot be zero.");
    }

    // --- PHASE 1: TOKENIZATION ---
    std::vector<int64_t> aa_indices;
    aa_indices.reserve(seq_length);

    for (char aa : peptide_sequence) {
        auto it = aa_vocab.find(std::toupper(aa)); // Ensure uppercase safety
        if (it == aa_vocab.end()) {
            throw std::invalid_argument(std::string("Unknown amino acid residue encountered: ") + aa);
        }
        aa_indices.push_back(it->second);
    }

    // 1. Sequence Tensor (using our newly mapped aa_indices)
    std::vector<int64_t> aa_shape = {batch_size, seq_length};
    Ort::Value aa_tensor = Ort::Value::CreateTensor<int64_t>(
        memory_info_, const_cast<int64_t*>(aa_indices.data()), aa_indices.size(), aa_shape.data(), aa_shape.size());

    // 2. Modification Tensor (109-Dimension)
    int64_t mod_dim = 109;
    std::vector<float> mod_x_data(batch_size * seq_length * mod_dim, 0.0f);
    std::vector<int64_t> mod_x_shape = {batch_size, seq_length, mod_dim};
    Ort::Value mod_tensor = Ort::Value::CreateTensor<float>(
        memory_info_, mod_x_data.data(), mod_x_data.size(), mod_x_shape.data(), mod_x_shape.size());

    // 3. Charge Tensor (2D: [batch_size, 1])
    std::vector<float> charge_data = {charge};
    std::vector<int64_t> charge_shape = {batch_size, 1};
    Ort::Value charge_tensor = Ort::Value::CreateTensor<float>(
        memory_info_, charge_data.data(), charge_data.size(), charge_shape.data(), charge_shape.size());

    // 4. NCE Tensor (2D: [batch_size, 1])
    std::vector<float> nce_data = {nce};
    std::vector<int64_t> nce_shape = {batch_size, 1};
    Ort::Value nce_tensor = Ort::Value::CreateTensor<float>(
        memory_info_, nce_data.data(), nce_data.size(), nce_shape.data(), nce_shape.size());

    // 5. Instrument Index Tensor (1D: [batch_size])
    std::vector<int64_t> inst_data = {instrument_index};
    std::vector<int64_t> inst_shape = {batch_size};
    Ort::Value inst_tensor = Ort::Value::CreateTensor<int64_t>(
        memory_info_, inst_data.data(), inst_data.size(), inst_shape.data(), inst_shape.size());

    // Execute exactly as before...
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

} // namespace OpenMS