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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

#include <torch/torch.h>
#include <torch/script.h>
#include <vector>
#include <iostream>
#include <OpenMS/ML/AlphapeptdeepWrapper.h>
#include <OpenMS/ML/AlphaDatahandling.h>
using namespace OpenMS;

#include <test_config.h>

START_TEST(Torch_Integration, "$Id$")

//for Testing purposes

std::string model_path =OPENMS_ML_GET_TEST_DATA_PATH("models/serialized_model_script.zip");
std::string model_config_path =OPENMS_ML_GET_TEST_DATA_PATH("/models/model.model_const.txt");

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


END_TEST