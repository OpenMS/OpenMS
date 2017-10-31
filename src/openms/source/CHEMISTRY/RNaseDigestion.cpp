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
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>

using namespace std;

namespace OpenMS
{
  void RNaseDigestion::setEnzyme(const String& enzyme_name)
  {
    enzyme_ = RNaseDB::getInstance()->getEnzyme(enzyme_name);
  }

  void RNaseDigestion::digest(const String& rna, vector<String>& output,
                              Size min_length, Size max_length) const
  {
    // initialization
    output.clear();
    if (rna.empty() || ((rna.size() == 1) && (rna[0] == 'p'))) return;

    // handle terminal phosphates in input:
    bool has_5prime_p = (rna[0] == 'p');
    bool has_3prime_p = (rna[rna.size() - 1] == 'p');
    // mark sequence ends:
    String temp_rna = "^" + rna.substr(has_5prime_p, rna.size() - has_5prime_p -
                                       has_3prime_p) + "$";

    // beginning positions of "naive" fragments:
    vector<Size> fragment_pos = tokenize_(temp_rna);
    // after "^" or before "$" aren't valid cleavages:
    if (fragment_pos.size() > 1)
    {
      if (fragment_pos[1] == 1)
      {
        fragment_pos.erase(++fragment_pos.begin());
      }
      if (fragment_pos[fragment_pos.size() - 1] == temp_rna.size() - 1)
      {
        fragment_pos.resize(fragment_pos.size() - 1);
      }
    }

    vector<StringView> unmod_output;
    // don't apply length filters yet, because we modified the original string:
    digestAfterTokenize_(fragment_pos, temp_rna, unmod_output);

    String three_prime_gain =
      dynamic_cast<const DigestionEnzymeRNA*>(enzyme_)->getThreePrimeGain();
    String five_prime_gain =
      dynamic_cast<const DigestionEnzymeRNA*>(enzyme_)->getFivePrimeGain();

    for (vector<StringView>::iterator it = unmod_output.begin();
         it != unmod_output.end(); ++it)
    {
      String fragment = it->getString();
      Size actual_length = fragment.size();
      bool is_5prime_end = (fragment[0] == '^');
      bool is_3prime_end = (fragment[fragment.size() - 1] == '$');
      if (is_5prime_end) // original 5' end -> no 5' enzyme mod
      {
        actual_length--; // don't count the "^"
        if (has_5prime_p)
        {
          fragment[0] = 'p';
        }
        else
        {
          fragment = fragment.substr(1);
        }
      }
      else
      {
        fragment = five_prime_gain + fragment;
      }
      if (is_3prime_end) // original 3' end -> no 3' enzyme mod
      {
        actual_length--; // don't count the "$"
        if (has_3prime_p)
        {
          fragment[fragment.size() - 1] = 'p';
        }
        else
        {
          fragment = fragment.substr(0, fragment.size() - 1);
        }
      }
      else
      {
        fragment += three_prime_gain;
      }

      if ((actual_length >= min_length) && (actual_length <= max_length))
      {
        output.push_back(fragment);
      }
    }
  }

} //namespace
