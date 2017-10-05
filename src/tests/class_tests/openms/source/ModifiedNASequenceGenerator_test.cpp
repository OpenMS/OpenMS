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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ModifiedNASequenceGenerator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModifiedNASequenceGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


START_SECTION((
	static void applyFixedModifications(
		const std::vector<ModifiedNASequenceGenerator::ConstRibonucleotidePtr>::const_iterator& fixed_mods_begin,
		const std::vector<ModifiedNASequenceGenerator::ConstRibonucleotidePtr>::const_iterator& fixed_mods_end,
		NASequence& sequence)))
{
  vector<ModifiedNASequenceGenerator::ConstRibonucleotidePtr> fixed_mods;
  NASequence sequence;

  ModifiedNASequenceGenerator::applyFixedModifications(
  	fixed_mods.begin(),
	  fixed_mods.end(),
	  sequence);
}
END_SECTION

START_SECTION(applyVariableModifications())
  vector<ModifiedNASequenceGenerator::ConstRibonucleotidePtr> var_mods;
  NASequence sequence;
  vector<NASequence> all_modified_sequences;

  ModifiedNASequenceGenerator::applyVariableModifications(
		var_mods.begin(),
  	var_mods.end(),
	  sequence,
	  1,
    all_modified_sequences,
	  true);

END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



