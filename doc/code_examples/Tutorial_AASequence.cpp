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

//! [AASequence]

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // generate AASequence object from String
  const String s = "DEFIANGER";
  AASequence peptide1 = AASequence::fromString(s);

  // generate AASequence object from string literal
  AASequence peptide2 = AASequence::fromString("PEPTIDER");

  // extract prefix and suffix
  AASequence prefix(peptide1.getPrefix(2));
  AASequence suffix(peptide1.getSuffix(3));
  cout << peptide1.toString() << " "
       << prefix << " "
       << suffix << endl;
  
  // create chemically modified peptide
  AASequence peptide_meth_ox = AASequence::fromString("PEPTIDESEKUEM(Oxidation)CER");
  cout << peptide_meth_ox.toString() << " "
       << peptide_meth_ox.toUnmodifiedString()
       << endl;

  // mass of the full, uncharged peptide
  double peptide_mass_mono = peptide_meth_ox.getMonoWeight();
  cout << "Monoisotopic mass of the uncharged, full peptide: " << peptide_mass_mono << endl;

  double peptide_mass_avg = peptide_meth_ox.getAverageWeight();
  cout << "Average mass of the uncharged, full peptide: " << peptide_mass_avg << endl;

  // mass of the 2+ charged b-ion with the given sequence
  double ion_mass_2plus = peptide_meth_ox.getMonoWeight(Residue::BIon, 2);
  cout << "Mass of the doubly positively charged b-ion: " << ion_mass_2plus << endl;
         
  // ... many more
  return 0;
}

//! [AASequence]
