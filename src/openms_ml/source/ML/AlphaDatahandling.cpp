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

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <OpenMS/ML/AlphaDatahandling.h>
#include <unordered_map>
using namespace OpenMS;


AlphaDatahandler::AlphaDatahandler(const std::vector<std::string>& seq)
{
    sequences = seq;
    
    //import dummy data
    std::vector<int> temp_mod_sites;
    temp_mod_sites.push_back (3);
    temp_mod_sites.push_back (5);
    
    std::string ox = "Oxidation@M";
    std::string car = "Carbodimethylation@C";
    
    std::vector<std::string> temp_mods;
    temp_mods.push_back (ox);
    temp_mods.push_back (car);
    
    for (int i = 0; i < sequences.size(); ++i)
    {
        const std::string& seq = sequences[i];
        nAA.push_back(seq.length());
        retention_time.push_back(727.0444);
        mod_sites.push_back(temp_mod_sites);
        mods.push_back(temp_mods);
    }

}

void AlphaDatahandler::printAll()
{
    for (int i = 0; i < sequences.size(); ++i)
    {
        std::cout<<"sequence: "<<sequences[i]<<" len: "<<nAA[i]<<" RT: "<<retention_time[i]<<"\n";
        
        std::cout<<"mod_sites: ";
        for (const auto& row : mod_sites[i]) 
        {
                std::cout << row << " ";
        }
        
        std::cout<<"\nmod_sites: ";
        for (const auto& row : mods[i]) 
        {
                std::cout << row << " ";
            
        }
        std::cout<<"\n";
    }
    
}


std::vector<std::vector<long int>> AlphaDatahandler::get_batch_aa_indices()
{
    // Convert peptide sequences into AA ID array. ID=0 is reserved for masking,
    // so ID of 'A' is 1, ID of 'B' is 2, ..., ID of 'Z' is 26 (maximum).
    // Zeros are padded for the N- and C-term modifications for each sequence.

    std::vector<std::vector<long int>> aa_indices(sequences.size(), std::vector<long int>(sequences[0].size() + 2));

    for (int i = 0; i < sequences.size(); ++i)
    {
        const std::string& seq = sequences[i];
        std::vector<long int>& indices = aa_indices[i];

        indices[0] = 0;
        for (int j = 0; j < seq.size(); ++j)
        {
            indices[j + 1] = (long int)seq[j] - 'A' + 1;
        }
        indices[seq.size() + 1] = 0;
    }

    return aa_indices;
}

std::vector<std::vector<std::vector<int>>> AlphaDatahandler::get_batch_mod_feature(const std::vector<std::string>& mod_elements)
{
    int len_mod_elements = mod_elements.size();
    int seq_len = nAA[0]+2;

    std::vector<std::vector<std::vector<int>>> mod_x_batch(sequences.size(), std::vector<std::vector<int>>(seq_len, std::vector<int>(len_mod_elements, 0.0)));
    
    std::unordered_map<std::string, int> mod_elem_to_idx;
    for (int i=0; i<mod_elements.size(); ++i) 
    {
        mod_elem_to_idx[mod_elements[i]] = 0;
    }
    
    //dummy CHNOPS values
    std::vector<int> my_vec( mod_elements.size(), 0.0);
    my_vec[0] = 11.0;
    my_vec[1] = 3.0;
    
    for (int i=0; i<mod_x_batch.size();++i)
    {
      if(mod_sites[i].size()>0)
      {
        for (int site: mod_sites[i]) 
        {
            mod_x_batch[i][site] = my_vec;
        }
      }
    }
    
    /*for (int i=0; i<mod_x_batch.size(); ++i)
    {
      for(int j=0; j<mod_x_batch[i].size(); ++j)
      {
        std::cout<<i<<","<<" [";
        for (int k=0; k<mod_x_batch[i][j].size(); ++k) 
        {
            std::cout<<mod_x_batch[i][j][k]<<" ";
        }
        std::cout<<"]\n";
      }
    }*/
    
    
    return mod_x_batch;
}
