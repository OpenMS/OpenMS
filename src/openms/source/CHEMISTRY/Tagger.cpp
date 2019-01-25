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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/Tagger.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

namespace OpenMS
{
  char Tagger::getAAByMass_(float m) const
  {
     // fast check for border cases
     if (m < min_gap_ || m > max_gap_) return ' ';

     auto left = mass2aa.lower_bound(m - m * ppm_ * 1e-6);
     auto right = mass2aa.lower_bound(m + m * ppm_ * 1e-6);
     if (left == right) return ' '; // not found in tolerance window
     return left->second; // select first match
  }               

  void Tagger::getTag_(std::string & tag, const std::vector<float>& mzs, const size_t i, std::set<std::string>& tags) const 
  {
    const size_t N = mzs.size();
    size_t j = i + 1;
    // recurse for all peaks in distance < max_gap
    while (j < N) 
    {
      const float gap = mzs[j] - mzs[i];
      //cout << i << "\t" << j << "\t" << mzs[i] << "\t" << mzs[j] << "\t" << gap << endl;
      if (gap > max_gap_) { return; } // already too far away - continue with next parent
      if (tag.size() == max_tag_length_) { return; } // maximum tag size reached? - continue with next parent
      const char aa = getAAByMass_(gap);
      if (aa == ' ') { ++j; continue; } // can't extend tag
      tag += aa;
      getTag_(tag, mzs, j, tags);
      if (tag.size() < min_tag_length_) { tag.resize(tag.size() - 1); ++j; continue; }
      tags.insert(tag);
      tag.resize(tag.size() - 1);  // remove last char
      ++j;
    }         
  }

  Tagger::Tagger(size_t min_tag_length, float ppm, size_t max_tag_length)
  {
    ppm_ = ppm;
    min_tag_length_ = min_tag_length;
    max_tag_length_ = max_tag_length;
    const std::set<const Residue*> aas = ResidueDB::getInstance()->getResidues();
    for (auto r : aas)
    {
      const char letter = r->getOneLetterCode()[0]; 
      const float mass = r->getMonoWeight(Residue::Internal);
      mass2aa[mass] = letter;
      //cout << "Mass: " << mass << "\t" << letter << endl; 
    }
    min_gap_ = mass2aa.begin()->first - mass2aa.begin()->first * ppm * 1e-6;
    max_gap_ = mass2aa.rbegin()->first + mass2aa.rbegin()->first * ppm * 1e-6;
  }

  void Tagger::getTag(const std::vector<float>& mzs, std::set<std::string>& tags) const 
  {
    // start peak
    std::string tag;
    tag.reserve(30);
    for (size_t i = 0; i < mzs.size() - min_tag_length_; ++i)
    {
      getTag_(tag, mzs, i, tags);
      tag.clear();
    }
  }

  void Tagger::getTag(const MSSpectrum& spec, std::set<std::string>& tags) const
  {
    const size_t N = spec.size();
    if (N < min_tag_length_) { return; }
    // copy to float vector (speed)
    std::vector<float> mzs;
    mzs.reserve(N);
    for (auto const & p : spec) { mzs.push_back(p.getMZ()); }
    getTag(mzs, tags); 
  }
}

