//! [ResidueModification]
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

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // construct a AASequence object, query a residue
  // and output some of its properties
  AASequence aas = AASequence::fromString("DECIANGER");
  cout << aas[2].getName() << " "
       << aas[2].getFormula().toString() << " "
       << aas[2].getModificationName() << " "
       << aas[2].getMonoWeight() << endl;
  
  // find a modification in ModificationsDB
  // and output some of its properties
  // getInstance() returns a pointer to a ModsDB instance
  ResidueModification mod = ModificationsDB::getInstance()->getModification("Carbamidomethyl (C)");
  cout << mod.getOrigin() << " "
       << mod.getFullId() << " "
       << mod.getDiffMonoMass() << " "
       << mod.getMonoMass() << endl;
  
  // set the modification on a residue of a peptide
  // and output some of its properties (the formula and mass have changed)
  // in this case ModificationsDB is used in the background
  // to relate the name of the mod to its attributes
  aas.setModification(2, "Carbamidomethyl (C)");
  cout << aas[2].getName() << " "
   	<< aas[2].getFormula().toString() << " "
   	<< aas[2].getModificationName() << " "
   	<< aas[2].getMonoWeight() << endl;

  return 0;
} //end of main

//! [ResidueModification]
