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
    FASTAFile::load(fasta, fastaEntry);
    for (auto &fe : fastaEntry)
    {
      auto aaSeq = AASequence::fromString(fe.sequence);
      TheoreticalSpectrumGenerator generator;
      PeakSpectrum spec;
      generator.getSpectrum(spec, aaSeq, 1, 1);
      proteinVectors.push_back(spec);
    }
    #pragma omp parallel num_threads(3)
    {
      std::cout << "Total " << proteinVectors.size() << " proteins processed\n";
    }
  }

  FLASHProFilterAlgorithm::~FLASHProFilterAlgorithm()
  {

  }//-march=native

  /*
      Size m = *specVectorIndex.rbegin();
      evergreen::Tensor<evergreen::cpx> y({m+1});
      for (auto j: specVectorIndex)
        y[j] = evergreen::cpx{double(specVector[j]), .0};
      for(int i=0;i<proteinVectors.size();i++)
      {
        Size  n = *proteinVectorIndex[i].rbegin();
        evergreen::Tensor<evergreen::cpx> x({n+1});

        //std::cout<<n+1<<std::endl;

        for (auto &j: proteinVectorIndex[i])
        {
         // std::cout << j << " "<< i<< std::endl;
          x[j] = evergreen::cpx{double(proteinVectors[i][j]), .0};
        }
        evergreen::Tensor<evergreen::cpx> z = fft_convolve(x, y);
      }

  */

  std::map<int, double> FLASHProFilterAlgorithm::getScores(MSSpectrum &decovSpec, double intThreshold)
  {

    MSSpectrum filtered;
    filtered.reserve(decovSpec.size());
    for (auto &sp : decovSpec)
    {
      if (sp.getIntensity() <= intThreshold)
      {
        continue;
      }
      filtered.push_back(Peak1D(sp.getMZ(), log10(sp.getIntensity())));
    }

    std::map<int, double> scores;
    if (filtered.size() == 0)
    {
      return scores;
    }

    auto maxPeakMass = filtered[filtered.size() - 1].getMZ();
    auto size = FLASHDeconvAlgorithm::getNominalMass(maxPeakMass) + 1;


    //std::cout<<"threads="<<omp_get_max_threads()<<std::endl;


    #pragma omp parallel for
    for (int i = 0; i < proteinVectors.size(); i++)
    {
      auto &pSpec = proteinVectors[i];
      auto maxMz = pSpec[pSpec.size() - 1].getMZ();
      auto convSize = size + FLASHDeconvAlgorithm::getNominalMass(maxMz) + 2;
      auto *vector = new float[convSize];
      std::fill_n(vector, convSize, 0);
      boost::dynamic_bitset<> indices = boost::dynamic_bitset<>(convSize);

      for (auto &pp : pSpec)
      {
        for (auto &sp : filtered)
        {
          auto sumMz = pp.getMZ() + sp.getMZ();
          auto loc = FLASHDeconvAlgorithm::getNominalMass(sumMz);
          vector[loc] += pp.getIntensity() * sp.getIntensity();
          indices[loc] = true;
        }
      }
      delete[] vector;
    }
    return scores;
  }


}