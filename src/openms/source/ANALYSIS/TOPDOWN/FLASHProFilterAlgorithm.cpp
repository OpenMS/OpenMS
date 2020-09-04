// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHProFilterAlgorithm.h>

namespace OpenMS
{
  FLASHProFilterAlgorithm::FLASHProFilterAlgorithm(const String& fasta)
  {
    FASTAFile::load(fasta,fastaEntry);
    for (auto& fe : fastaEntry)
    {
      auto aaSeq = AASequence::fromString(fe.sequence);
      TheoreticalSpectrumGenerator generator;
      PeakSpectrum spec;
      generator.getSpectrum(spec, aaSeq, 1, 1);
      Byte* pvector = nullptr;
      std::set<Size> pindices;
      specToVectors(spec, pindices, pvector);
      proteinVectors.push_back(pvector);
      proteinVectorIndex.push_back(pindices);
    }
    std::cout<<"Total "<<proteinVectors.size()<< " proteins\n";
  }

  FLASHProFilterAlgorithm::~FLASHProFilterAlgorithm(){
    for(auto & pvector : proteinVectors){
      delete[] pvector;
    }
  }

  std::map<int, double> FLASHProFilterAlgorithm::getScores(MSSpectrum &decovSpec)
  {
    std::set<Size> specVectorIndex;
    Byte *specVector = nullptr;
    specToVectors(decovSpec, specVectorIndex, specVector);

    std::cout<<*specVectorIndex.rbegin()<<std::endl;//TODO

    std::map<int, double> scores;
    for(int i=0;i<proteinVectors.size();i++){
      auto size = *proteinVectorIndex[i].rbegin() + *specVectorIndex.rbegin(); // check
      Byte* vector = new Byte[size];
      std::fill_n(vector, size, 0);
      std::set<Size> indices;
      for (auto &pi : proteinVectorIndex[i])
      {
        auto pv = proteinVectors[i][pi];
        for (auto &si : specVectorIndex)
        {
          auto sv = specVector[si];
          vector[pi + si] += pv * sv;
          indices.insert(pi + si);
        }
      }
    }

    delete[] specVector;

    return scores;
  }

  void FLASHProFilterAlgorithm::specToVectors(MSSpectrum &spec,
                                              std::set<Size> &vindex,
                                              Byte * vector)
  {
    auto maxPeakMass = spec[spec.size()-1].getMZ();
    auto size = FLASHDeconvAlgorithm::getNominalMass(maxPeakMass);
    std::fill_n(vector, size, 0);
    for(auto & p : spec){
      auto pm = p.getMZ();
      Size index = FLASHDeconvAlgorithm::getNominalMass(pm);
      vector[index] ++;
      vindex.insert(index);
    }

  }
}