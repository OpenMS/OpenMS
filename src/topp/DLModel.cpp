#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <onnxruntime_cxx_api.h>

class DLModel
{
public:
  int run(const std::string& model_path,
          const std::string& csv_path,
          const std::string& output_path)
  {
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "DLModel");
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);
    Ort::Session session(env, model_path.c_str(), session_options);

    // ========= 1. Parse CSV =========
    std::ifstream file(csv_path);
    if (!file.is_open())
    {
      std::cerr << "Could not open CSV: " << csv_path << std::endl;
      return 1;
    }

    std::string line;
    bool header = true;

    std::vector<int> feature_ids;
    std::vector<std::vector<std::vector<float>>> groups;

    int current_id = -1;
    std::vector<std::vector<float>> current_group;

    while (std::getline(file, line))
    {
      if (header) { header = false; continue; }
      std::stringstream ss(line);
      std::string cell;
      std::vector<std::string> tokens;
      while (std::getline(ss, cell, ','))
        tokens.push_back(cell);

      int feature_idx = std::stoi(tokens[0]);
      std::vector<float> feat;
      for (size_t i = 1; i < tokens.size() - 1; i++) // exclude FeatureIndex + Class
        feat.push_back(std::stof(tokens[i]));

      // start new group if FeatureIndex changes
      if (feature_idx != current_id && current_id != -1)
      {
        groups.push_back(current_group);
        current_group.clear();
      }
      if (feature_idx != current_id)
      {
        feature_ids.push_back(feature_idx);
        current_id = feature_idx;
      }
      current_group.push_back(feat);
    }
    if (!current_group.empty())
      groups.push_back(current_group);

    size_t feat_dim = groups[0][0].size();
    std::cout << "Loaded " << groups.size() << " FeatureIndex groups" << std::endl;

    Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(
        OrtDeviceAllocator, OrtMemTypeCPU);

    const char* input_names[] = {"input"};
    const char* output_names[] = {"output"};

    std::ofstream fout(output_path);
    fout << "FeatureIndex,Prediction\n";

    // ========= 2. Run inference per group =========
    for (size_t g = 0; g < groups.size(); g++)
    {
      auto& group = groups[g];
      size_t rows = group.size();

      // flatten group into contiguous buffer
      std::vector<float> flat_features;
      flat_features.reserve(rows * feat_dim);
      for (auto& row : group)
        flat_features.insert(flat_features.end(), row.begin(), row.end());

      // input tensor [rows, feat_dim]
      std::array<int64_t, 2> input_shape{(int64_t)rows, (int64_t)feat_dim};
      Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
          memory_info, flat_features.data(), flat_features.size(),
          input_shape.data(), input_shape.size());

      std::vector<Ort::Value> input_tensors;
      input_tensors.push_back(std::move(input_tensor));

      auto output_tensors = session.Run(
          Ort::RunOptions{nullptr},
          input_names,
          input_tensors.data(),
          input_tensors.size(),
          output_names,
          1);

      float* output_data =
          output_tensors.front().GetTensorMutableData<float>();

      // Model should already aggregate across rows → output is scalar
      fout << feature_ids[g] << "," << output_data[0] << "\n";
    }

    fout.close();
    return 0;
  }
};

// ========== Main entry ==========
int main(int argc, char* argv[])
{
  if (argc != 4)
  {
    std::cerr << "Usage: " << argv[0]
              << " <model.onnx> <input.csv> <output.csv>\n";
    return 1;
  }

  DLModel model;
  return model.run(argv[1], argv[2], argv[3]);
}
