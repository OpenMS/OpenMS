// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
//  --------------------------------------------------------------------------

#include <OpenMS/ML/AlphapeptdeepWrapper.h>
#include <OpenMS/ML/AlphaDatahandling.h>
#include <OpenMS/ML/OpenMS_MLConfig.h>
#include <torch/torch.h>
#include <torch/script.h>
#include <vector>
using namespace OpenMS;

struct Alphapeptdeepwrapper::Impl 
{
    torch::jit::script::Module model;
    config_param config_data;
  
    Impl(const std::string& model_path, const std::string& model_config_path)
    {
       try 
       {
        	model = torch::jit::load(model_path);
        }
        catch (const c10::Error& e) 
        {
          std::cerr << "Error loading the model..\n";
        }
        
        config_data = Alphapeptdeepwrapper::set_config_param(model_config_path);
        
    };
    
    ~Impl()=default;
  
};

 
Alphapeptdeepwrapper::Alphapeptdeepwrapper(const std::string& model_path, const std::string& model_config_path)
    : pimpl(new Impl(model_path, model_config_path)){}
 
Alphapeptdeepwrapper::~Alphapeptdeepwrapper() = default;
 
config_param Alphapeptdeepwrapper::set_config_param(const std::string& model_config_path)
{
    config_param data;
    std::ifstream infile(model_config_path);
    std::string line;
    
    std::string check_key = "mod_elements:";
    while (getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string key, value;
        ss >> key >> value;
        
        if (check_key == "mod_elements:")
        {
            if (key == "-")
            { 
               data.mod_elements.push_back(value);
            }
            else
            {
              check_key = key;
            }
        } 
        else if (check_key == "instruments:") 
        {
            if (key == "-")
            { 
               data.instruments.push_back(value);
            }
            else
            {
               check_key = key;
            }
            
        } 
        
        if (key == "max_instrument_num:") 
        {
            data.max_instrument_num = std::stoi(value);
        } 
        else if (key == "aa_embedding_size:") 
        {  
            data.aa_embedding_size = std::stoi(value);
        }
       
    }
    
    return data;
}


std::vector<float> Alphapeptdeepwrapper::predict(const std::vector<std::string>& seq)
{  
    AlphaDatahandler Data_handler(seq);
    
    // Generate aa_indices
    std::vector<std::vector<long int>> aa_indices = Data_handler.get_batch_aa_indices();
   
    // convert 2D vector (aa_indices) to torch tensor
    std::vector<long int> aa_indices_1d;
    for (int i = 0; i < (static_cast<long>(aa_indices.size())); i++) {
        aa_indices_1d.insert(aa_indices_1d.end(), aa_indices[i].begin(), aa_indices[i].end());
    }
  
    long int* aa_indices_data = aa_indices_1d.data();
  
    torch::Tensor aa_indices_tensor = torch::from_blob(
        aa_indices_data,                      
        {static_cast<long>(aa_indices.size()), static_cast<long>(aa_indices[0].size())},   
        torch::kLong                           
    );
      
    
    // Generate mod_matrix 
    std::vector<std::vector<std::vector<int>>> mod_x_batch =  Data_handler.get_batch_mod_feature(pimpl->config_data.mod_elements);   
   
    // convert 3D vector (mod_x_batch) to torch tensor
    std::vector<int> mod_x_batch_1d;
    for (int i = 0; i < (static_cast<long>(mod_x_batch.size())); i++) {
        for (int j = 0; j < (static_cast<long>(mod_x_batch[i].size())); j++) {
            mod_x_batch_1d.insert(mod_x_batch_1d.end(), mod_x_batch[i][j].begin(), mod_x_batch[i][j].end());
        }
    }
  
    int* mod_x_batch_data = mod_x_batch_1d.data();
  
    torch::Tensor mod_x_batch_tensor = torch::from_blob(
        mod_x_batch_data,                       
        {static_cast<long>(mod_x_batch.size()), static_cast<long>(mod_x_batch[0].size()), static_cast<long>(mod_x_batch[0][0].size())}, 
        torch::kInt                              
    );
  
    
    torch::Tensor mod_x_batch_tensor_f = mod_x_batch_tensor.toType(torch::kFloat);
    
    //create inputs as model accepted
    std::vector<torch::jit::IValue> features;
        features.push_back(aa_indices_tensor);
        features.push_back(mod_x_batch_tensor_f);
        
    //predict from model :already in eval mode
    torch::jit::IValue output = pimpl->model.forward(features);
    torch::Tensor output_tensor = output.toTensor();
    std::vector<float> output_vector(output_tensor.data_ptr<float>(), output_tensor.data_ptr<float>() + output_tensor.numel());
  
    return output_vector;
}






