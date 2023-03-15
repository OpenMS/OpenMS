#include <torch/torch.h>
#include <torch/script.h>
#include <vector>
#include <OpenMS/ML/AlphapeptdeepWrapper.h>
#include <OpenMS/ML/AlphaDatahandling.h>
using namespace OpenMS;

int main()
{
  //for Testing purposes
  
  std::string model_path = "../models/serialized_model_script.zip";
  std::string model_config_path = "../models/model.model_const.txt";
  
  Alphapeptdeepwrapper TestModel(model_path, model_config_path);
  
  //seq should be the same size and pass modified sequence so, can extract everything in constructor of handler
  std::vector<std::string> seq_array;
  seq_array.push_back("AAAALAGGKKSK");
  seq_array.push_back("EMMSQVTLQHMN");
    
  std::vector<float> rt_pred = TestModel.predict(seq_array);
  
  std::cout<<"Predictions.. \n";
  for (float res: rt_pred) 
  {
      std::cout<<res<<" ";
  }
  
  
}
