// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Alexandra Scherbart $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/PIP/PeakIntensityPredictor.h>

using namespace OpenMS;
using namespace std;

Int main()
{
  //Create a vector for the predicted values that is large enough to hold them all
  vector<AASequence> peptides;
  peptides.push_back(AASequence("IVGLMPHPEHAVEK"));
  peptides.push_back(AASequence("LADNISNAMQGISEATEPR"));
  peptides.push_back(AASequence("ELDHSDTIEVIVNPEDIDYDAASEQAR"));
  peptides.push_back(AASequence("AVDTVR"));
  peptides.push_back(AASequence("AAWQVK"));
  peptides.push_back(AASequence("FLGTQGR"));
  peptides.push_back(AASequence("NYPSDWSDVDTK"));
  peptides.push_back(AASequence("GSPSFGPESISTETWSAEPYGR"));
  peptides.push_back(AASequence("TELGFDPEAHFAIDDEVIAHTR"));

  //Create new predictor model with vector of AASequences
  PeakIntensityPredictor model;

  //Perform prediction with LLM model
  vector<DoubleReal> predicted = model.predict(peptides);

  //for each element in peptides print sequence as well as corresponding predicted peak intensity value.
  for (Size i = 0; i < peptides.size(); i++)
  {
    cout << "Intensity of " << peptides[i] << " is " << predicted[i] << endl;
  }

  return 0;
} //end of main
