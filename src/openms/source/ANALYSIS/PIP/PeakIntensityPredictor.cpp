// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/PIP/PeakIntensityPredictor.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <cmath>
#include <string>

using namespace std;


namespace OpenMS
{

  PeakIntensityPredictor::PeakIntensityPredictor() :
    llm_()
  {
  }

  PeakIntensityPredictor::~PeakIntensityPredictor()
  {
  }

  double PeakIntensityPredictor::predict(const AASequence & sequence)
  {
    //get corresponding attribute values
    vector<double> aafeat = getPropertyVector_(sequence);
    //normalize vector (center and scale by variance)
    llm_.normalizeVector(aafeat);
    //pass data_ to llm model and return predicted value
    return map_(aafeat);
  }

  double PeakIntensityPredictor::predict(const AASequence & sequence, vector<double> & add_info)
  {
    //get corresponding attribute values
    vector<double> aafeat = getPropertyVector_(sequence);
    //normalize vector (center and scale by variance)
    llm_.normalizeVector(aafeat);
    //pass data_ to llm model and return predicted value
    double tmp = map_(aafeat);

    //Determine corresponding winning prototype and number of cluster
    //peptide is assigned to. Additionally, the error is given.
    add_info = calculateAddInfo_(aafeat);

    return tmp;
  }

  vector<double> PeakIntensityPredictor::predict(const vector<AASequence> & sequences)
  {
    vector<double> out(sequences.size());

    for (Size i = 0; i < sequences.size(); i++)
    {
      out[i] = predict(sequences[i]);
    }

    return out;
  }

  vector<double> PeakIntensityPredictor::predict(const vector<AASequence> & sequences, vector<vector<double> > & add_info)
  {
    vector<double> out(sequences.size());
    add_info.resize(sequences.size());

    //for each element in sequences
    for (Size i = 0; i < sequences.size(); i++)
    {
      out[i] = predict(sequences[i], add_info[i]);
    }

    return out;
  }

  double PeakIntensityPredictor::map_(const vector<double> & data)
  {
    double sum_g_i, g_i;
    Matrix<double> code = llm_.getCodebooks();
    vector<double> wout = llm_.getVectorWout();
    Matrix<double> A = llm_.getMatrixA();

    //determine best matching unit
    Size winner = findWinner_(data);
    //calculate Gaussian neighborhood function
    vector<double> nei = llm_.neigh(llm_.getCord(), winner, llm_.getLLMParam().radius);

    sum_g_i = 0.0;
    //for each prototype
    for (Size c = 0; c < code.rows(); c++)
    {
      //sum up total sum of nei
      sum_g_i += nei[c];
    }
    g_i = 0.0;

    //map data to best matching unit of prototypes and get
    //linear mapping
    for (Size i = 0; i < code.rows(); i++)
    {
      double c_x = 0.0;
      for (Size c = 0; c < code.cols(); c++)
      {
        c_x += (A.getValue(i, c) * (data[c] - code.getValue(i, c)));
      }
      //add linear bias to mapped value
      g_i += (c_x + wout[i]) * nei[i];
    }
    //divide by total neighborhood values
    g_i = g_i / sum_g_i;

    //normalize predicted values to distribution of training data
    g_i -= 3.364288;     //mean(targets)
    g_i /= 1.332298;     //sd(targets)

    return g_i;
  }

  Size PeakIntensityPredictor::findWinner_(const vector<double> & data)
  {
    Size winner = 0;
    double min_dist = 0.0;
    Matrix<double> code = llm_.getCodebooks();

    //calculate euclidean distance of vector data to prototype no 0.
    for (Size c = 0; c < data.size(); c++)
    {
      min_dist += (data[c] - code.getValue(0, c)) * (data[c] - code.getValue(0, c));
    }

    //calculate euclidean distance of vector data to the remaining prototypes
    for (Size i = 1; i < code.rows(); i++)
    {
      double dd = 0.0;
      for (Size c = 0; c < data.size(); c++)
      {
        dd += (data[c] - code.getValue(i, c)) * (data[c] - code.getValue(i, c));
      }
      //which is min?
      if (dd < min_dist)
      {
        winner = i;
        min_dist = dd;
      }
    }

    //return number of winning prototype with minimum distance to data
    return winner;
  }

  vector<double> PeakIntensityPredictor::calculateAddInfo_(const vector<double> & data)
  {
    vector<double> foo(3);
    Size winner = findWinner_(data);
    Matrix<double> code = llm_.getCodebooks();
    Matrix<UInt> cord = llm_.getCord();

    foo[0] = cord.getValue(winner, 0);
    foo[1] = cord.getValue(winner, 1);

    double dd = 0.0;
    //get distance of data to best matching unit (winner).
    for (Size c = 0; c < data.size(); c++)
    {
      dd += (data[c] - code.getValue(winner, c)) * (data[c] - code.getValue(winner, c));
    }
    //store the minimum distance in third column
    foo[2] = sqrt(dd);

    return foo;
  }

  std::vector<double> PeakIntensityPredictor::getPropertyVector_(const AASequence & sequence)
  {
    std::vector<double> out(18);

    //for each element in sequence = residue
    for (Size pos = 0; pos < sequence.size(); pos++)
    {
      char amino = sequence[pos].getOneLetterCode()[0];

      // numResidues of R
      out[0] += (amino == 'R' ? 1.0 : 0.0);
      //The Kerr-constant increments
      out[7] += AAIndex::getKHAG800101(amino);
      //Relative population of conformational state E
      out[15] += AAIndex::getVASM830103(amino);
      //Hydropathy scale (36% accessibility)
      out[10] += AAIndex::getNADH010106(amino);
      //Hydropathy scale (50% accessibility)
      out[11] += AAIndex::getNADH010107(amino);
      //Hydrophobicity coefficient in RP-HPLC, C8 with 0.1%TFA/MeCN/H2 O,
      out[16] += AAIndex::getWILM950102(amino);
      //Information measure for extended without H-bond,
      out[14] += AAIndex::getROBB760107(amino);
      //Optimized average non-bonded energy per atom,
      out[12] += AAIndex::getOOBM850104(amino);
      //Positive charge
      out[3] += AAIndex::getFAUJ880111(amino);
      //Helix-coil equilibrium constant
      out[4] += AAIndex::getFINA770101(amino);
      //Signal sequence helical potential
      out[1] += AAIndex::getARGP820102(amino);
      out[8] += (amino == 'M' ? 1.0 : 0.0);    // numResidues of M
      out[2] += (amino == 'F' ? 1.0 : 0.0);    // numResidues of F
      out[6] += (amino == 'H' ? 1.0 : 0.0);    // numResidues of H
      out[13] += (amino == 'Q' ? 1.0 : 0.0);    // numResidues of Q
      out[17] += (amino == 'Y' ? 1.0 : 0.0);    // numResiduesof Y
    }

    out[5] = AAIndex::calculateGB(sequence, 500.0);     //Estimated gas-phase basicity at 500 K
    out[9] = sequence.getAverageWeight();

    return out;
  }

}
