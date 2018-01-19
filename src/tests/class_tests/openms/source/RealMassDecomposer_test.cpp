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
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

using namespace OpenMS;
using namespace ims;
using namespace std;

Weights createWeights()
{
  Map<char, double> aa_to_weight;

  set<const Residue*> residues = ResidueDB::getInstance()->getResidues("Natural19WithoutI");

  for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
  {
    aa_to_weight[(*it)->getOneLetterCode()[0]] = (*it)->getMonoWeight(Residue::Internal);
  }

  // init mass decomposer
  IMSAlphabet alphabet;
  for (Map<char, double>::ConstIterator it = aa_to_weight.begin(); it != aa_to_weight.end(); ++it)
  {
    alphabet.push_back(String(it->first), it->second);
  }

  // initializes weights
  Weights weights(alphabet.getMasses(), 0.01);

  // optimize alphabet by dividing by gcd
  weights.divideByGCD();

  return weights;
}

START_TEST(RealMassDecomposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RealMassDecomposer* ptr = nullptr;
RealMassDecomposer* null_ptr = nullptr;

START_SECTION((RealMassDecomposer(const Weights &weights)))
{
  ptr = new RealMassDecomposer(createWeights());
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION


START_SECTION(~RealMassDecomposer())
{
	delete ptr;
}
END_SECTION


START_SECTION((decompositions_type getDecompositions(double mass, double error)))
{
  // TODO
}
END_SECTION

START_SECTION((decompositions_type getDecompositions(double mass, double error, const constraints_type &constraints)))
{
  // TODO
}
END_SECTION

START_SECTION((number_of_decompositions_type getNumberOfDecompositions(double mass, double error)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



