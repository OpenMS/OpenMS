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
  FLASHProFilterAlgorithm::FLASHProFilterAlgorithm(const String &fasta)
  {
    FASTAFile::load(fasta, fastaEntry);
    //std::set<int> map;
    for (auto &fe : fastaEntry)
    {
      auto aaSeq = AASequence::fromString(fe.sequence);
      TheoreticalSpectrumGenerator generator;
      PeakSpectrum spec;
      generator.getSpectrum(spec, aaSeq, 1, 1);
      proteinVectors.push_back(spec);
    }

    std::cout << "Total " << proteinVectors.size() << " proteins processed\n";
    //std::cout << map.size() << " " << *map.rbegin() << std::endl;
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

  std::vector<double> FLASHProFilterAlgorithm::getScores(MSSpectrum &decovSpec, double intThreshold)
  {
    int peakcntr = 150; //100 : 64, 50 : 38, 200 : 73
    static int cntr = 0;
    int window = 500;
    int count = 3;
    MSSpectrum filtered;
    filtered.reserve(peakcntr);
    std::vector<double> intensities;
    intensities.reserve(decovSpec.size());
    for (auto &sp : decovSpec)
    {
      if (sp.getIntensity() <= intThreshold)
      {
        continue;
      }
      intensities.push_back(sp.getIntensity());
      // filtered.push_back(Peak1D(sp.getMZ(), log10(sp.getIntensity())));
    }
    std::vector<double> scores;
    if (intensities.size() == 0)
    {
      return scores;
    }
    std::sort(intensities.rbegin(), intensities.rend());
    int th = intensities.size() > peakcntr ? peakcntr : intensities.size();
    double intthreshold2 = intensities[th];
    for (auto &sp : decovSpec)
    {
      if (sp.getIntensity() <= intthreshold2)
      {
        continue;
      }
      filtered.push_back(Peak1D(sp.getMZ(), log10(1 + sp.getIntensity()))); // log10
    }
    if (filtered.size() == 0 || decovSpec.getMSLevel() < 2)
    {
      return scores;
    }

    auto maxPeakMass = filtered[filtered.size() - 1].getMZ();
    auto size = FLASHDeconvAlgorithm::getNominalMass(maxPeakMass) + 1;
    //std::cout <<filtered.size()<<std::endl;

    scores.reserve(proteinVectors.size());
    #pragma omp parallel for
    for (int i = 0; i < proteinVectors.size(); i++)//
    {
      auto &pSpec = proteinVectors[i];

      auto maxMz = pSpec[pSpec.size() - 1].getMZ();
      auto convSize = size + FLASHDeconvAlgorithm::getNominalMass(maxMz) + 2;
      auto *vector = new float[convSize];
      std::fill_n(vector, convSize, 0);
      boost::dynamic_bitset<> indices = boost::dynamic_bitset<>(convSize);

      double max = .0;
      int maxIndex = 0;

      for (auto &pp : pSpec)
      {
        for (auto &sp : filtered)
        {
          auto sumMz = maxMz - pp.getMZ() + sp.getMZ();
          auto loc = FLASHDeconvAlgorithm::getNominalMass(sumMz);
          vector[loc] += pp.getIntensity() * sp.getIntensity();
          indices[loc] = true;
          if (max < vector[loc])
          {
            max = vector[loc];
            maxIndex = loc;
          }
        }
      }

      Size index = maxIndex - window;
      index = index < 0 ? 0 : index;
      std::vector<float> values;
      values.reserve(window * 2);
      while (index != indices.npos && index < maxIndex + window)
      {
        values.push_back(vector[index]);
        index = indices.find_next(index);
      }
      if (values.empty())
      {
        continue;
      }
      std::sort(values.rbegin(), values.rend());
      //std::cout <<max << ": " << values[0] << " " << values[1]<< " " << values[2] <<std::endl;
      int c = count < values.size() ? count : values.size();
      double score = std::accumulate(values.begin(), values.begin() + c, .0);

      scores[i] = score;
      delete[] vector;
    }
    if (!scores.empty())
    {
      double maxScore = *max_element(scores.begin(), scores.end());
      if (maxScore == scores[0])
      {
        cntr++;
      }
      std::cout << cntr << std::endl;
    }
    return scores;
  }
}