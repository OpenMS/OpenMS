// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <chrono>
#include <algorithm>

using namespace OpenMS;

DecoyGenerator::DecoyGenerator()
{
  const UInt64 seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  rng_ = boost::mt19937_64(seed);
}

void DecoyGenerator::setSeed(UInt64 seed)
{
  rng_.seed(seed);
}

AASequence DecoyGenerator::reverseProtein(const AASequence& protein) const
{
  OPENMS_PRECONDITION(!protein.isModified(), "Decoy generation only supports unmodified proteins.");
  String s = protein.toUnmodifiedString();
  std::reverse(s.begin(), s.end());
  return AASequence::fromString(s);
}

AASequence DecoyGenerator::reversePeptides(const AASequence& protein, const String& protease) const
{
  OPENMS_PRECONDITION(!protein.isModified(), "Decoy generation only supports unmodified proteins.");
  std::vector<AASequence> peptides;
  ProteaseDigestion ed;
  ed.setMissedCleavages(0); // important as we want to reverse between all cutting sites
  ed.setEnzyme(protease);
  ed.setSpecificity(EnzymaticDigestion::SPEC_NONE);
  ed.digest(protein, peptides);    
  String pseudo_reversed;
  for (int i = 0; i < static_cast<int>(peptides.size()) - 1; ++i)
  {
    std::string s = peptides[i].toUnmodifiedString();
    auto last = --s.end(); // don't reverse enzymatic cutting site
    std::reverse(s.begin(), last);
    pseudo_reversed += s;
  }
  // the last peptide of a protein is not an enzymatic cutting site so we do a full reverse
  std::string s = peptides[peptides.size() - 1 ].toUnmodifiedString();
  std::reverse(s.begin(), s.end());
  pseudo_reversed += s;
  return AASequence::fromString(pseudo_reversed);
}

AASequence DecoyGenerator::shufflePeptides(
        const AASequence& protein,
        const String& protease,
        const int max_attempts)
{
  OPENMS_PRECONDITION(!protein.isModified(), "Decoy generation only supports unmodified proteins.");

  std::vector<AASequence> peptides;
  ProteaseDigestion ed;
  ed.setMissedCleavages(0); // important as we want to reverse between all cutting sites
  ed.setEnzyme(protease);
  ed.setSpecificity(EnzymaticDigestion::SPEC_NONE);
  ed.digest(protein, peptides);    
  String protein_shuffled;
  for (int i = 0; i < static_cast<int>(peptides.size()) - 1; ++i)
  {
    std::string peptide_string = peptides[i].toUnmodifiedString();

    String peptide_string_shuffled = peptide_string;
    auto last = --peptide_string_shuffled.end();
    double lowest_identity(1.0);
    String lowest_identity_string(peptide_string_shuffled);
    for (int i = 0; i < max_attempts; ++i) // try to find sequence with low identity
    {
      shuffle_(std::begin(peptide_string_shuffled), last);

      double identity = SequenceIdentity_(peptide_string_shuffled, peptide_string);
      if (identity < lowest_identity)
      {
        lowest_identity = identity;
        lowest_identity_string = peptide_string_shuffled;

        if (identity <= (1.0/peptide_string_shuffled.size() + 1e-6)) break; // found perfect shuffle (only 1 (=cutting site) of all AAs match)
      }
    }
    protein_shuffled += lowest_identity_string;
  }
  // the last peptide of a protein is not an enzymatic cutting site so we do a full shuffle
  std::string peptide_string = peptides[peptides.size() - 1 ].toUnmodifiedString();
  String peptide_string_shuffled = peptide_string;
  double lowest_identity(1.0);
  String lowest_identity_string(peptide_string_shuffled);
  for (int i = 0; i < max_attempts; ++i) // try to find sequence with low identity
  {
    shuffle_(std::begin(peptide_string_shuffled), std::end(peptide_string_shuffled));
    double identity = SequenceIdentity_(peptide_string_shuffled, peptide_string);
    if (identity < lowest_identity)
    {
      lowest_identity = identity;
      lowest_identity_string = peptide_string_shuffled;
      if (identity == 0) break; // found best shuffle
    }
  }
  protein_shuffled += lowest_identity_string;
  return AASequence::fromString(protein_shuffled);
}

// static
double DecoyGenerator::SequenceIdentity_(const String& decoy, const String& target)
{
  int match = 0;
  for (Size i = 0; i < target.size(); ++i)
  {
    if (target[i] == decoy[i]) { ++match; }
  }
  double identity = (double) match / target.size();
  return identity;
}


